#!/usr/bin/env python3
"""Choose unused organisms, then reproducibly derive registered query workloads."""
import collections
from concurrent.futures import ThreadPoolExecutor

from scenario_common import *


def catalog():
    with (BASE / 'as.bacteria.txt').open() as f:
        for line in f:
            if line.startswith('#assembly_accession'):
                fields = line.lstrip('#').rstrip().split('\t'); break
        rows = list(csv.DictReader(f, fieldnames=fields, delimiter='\t'))
    lineage = {}
    for rank in ['genus', 'family', 'phylum']:
        lineage[rank] = {accession_root(r['assembly_accession']): r['genus_id']
                         for r in read_tsv(BASE / (rank + '_real_bacteria_speciesTaxid_origanismName'))}
    return rows, lineage


def previous_inputs(rows):
    used, manifests, paths = set(), {}, {}
    for pattern in ['*/input_sha256.json', '*/real_input_sha256.json']:
        for p in sorted((ROOT / 'runs').glob(pattern)):
            if p.parent == RUN: continue
            manifests[str(p)] = sha(p)
            for path in json.loads(p.read_text()):
                root = accession_root(path)
                if root:
                    used.add(root)
                    if Path(path).exists(): paths.setdefault(root, path)
    prior = read_tsv(OLD / 'natural_candidates.tsv')
    for r in prior:
        used.add(accession_root(r['object_id']))
        paths[accession_root(r['object_id'])] = r['path']
    # Paired accessions may have different numeric roots. Iterate to closure.
    while True:
        before = len(used)
        for r in rows:
            a, b = accession_root(r['assembly_accession']), accession_root(r['gbrs_paired_asm'])
            if a in used or b in used:
                used.add(a)
                if b: used.add(b)
        if len(used) == before: break
    samples = {r['biosample'] for r in rows if accession_root(r['assembly_accession']) in used
               and r['biosample'] not in ['', 'na']}
    cache = json.loads((OLD / 'prior_fingerprints.json').read_text())
    cache.update(json.loads((OLD / 'candidate_qc.json').read_text()))
    for r in prior: cache[r['path']] = r
    missing = [p for p in paths.values() if p not in cache]
    print('Additional prior sequence fingerprint checks:', len(missing), flush=True)
    with ThreadPoolExecutor(max_workers=32) as pool:
        for p, info in zip(missing, pool.map(inspect, missing)): cache[p] = info
    fingerprints = {v['fingerprint'] for v in cache.values()}
    return used, samples, fingerprints, manifests, {r['genus'] for r in prior}


def select_validation():
    complete = DATA / 'independence_audit.json'
    if complete.exists(): return read_tsv(DATA / 'validation_genomes.tsv')
    rows, lineage = catalog()
    used, samples, fingerprints, manifests, old_genera = previous_inputs(rows)
    files = {p.name.split('_', 2)[0] + '_' + p.name.split('_', 2)[1]: p
             for p in (BASE / 'bacteria').glob('GC[AF]*.fna.gz')}
    species, rejected = collections.defaultdict(list), collections.Counter()
    for r in rows:
        acc = r['assembly_accession']; root = accession_root(acc)
        if acc not in files: rejected['missing_local_version'] += 1; continue
        if root in used: rejected['previous_accession'] += 1; continue
        if r['biosample'] in samples or r['biosample'] in ['', 'na']:
            rejected['previous_or_missing_biosample'] += 1; continue
        if (r['group'] != 'bacteria' or r['genome_rep'] != 'Full' or
                r['excluded_from_refseq'] not in ['', 'na'] or
                not 1_000_000 <= int(r['genome_size']) <= 10_000_000 or int(r['contig_count']) > 300):
            rejected['metadata_qc'] += 1; continue
        if any(lineage[x].get(root, '') in ['0', '1', '2', 'na', ''] for x in lineage):
            rejected['missing_lineage'] += 1; continue
        if lineage['genus'][root] in old_genera: rejected['previous_genus'] += 1; continue
        row = dict(object_id=acc, path=str(files[acc]), species=r['species_taxid'],
                   biosample=r['biosample'], organism=r['organism_name'], assembly_level=r['assembly_level'],
                   **{x: lineage[x][root] for x in lineage})
        species[row['species']].append(row)
    genera = collections.defaultdict(list)
    for sid, values in species.items():
        by_sample = {}
        for r in sorted(values, key=lambda x: key(x['object_id'])): by_sample.setdefault(r['biosample'], r)
        species[sid] = list(by_sample.values())
        if len(by_sample) >= CONFIG['genomes_per_species']:
            genera[values[0]['genus']].append(sid)
    genera = {g: sorted(s, key=key) for g, s in genera.items() if len(s) >= CONFIG['species_per_genus']}
    print('Metadata-eligible unseen genera:', len(genera), flush=True)
    cachepath = DATA / 'candidate_qc.json'
    cache = json.loads(cachepath.read_text()) if cachepath.exists() else {}
    selected, chosen, selected_samples, selected_fp, failed = [], [], set(), set(), []
    phy, fam = collections.Counter(), collections.Counter()
    remaining = set(genera)
    while remaining and len(chosen) < CONFIG['validation_genera']:
        def order(g):
            r = species[genera[g][0]][0]
            return phy[r['phylum']], fam[r['family']], key(g)
        genus = min(remaining, key=order); remaining.remove(genus)
        group, group_species = [], []
        for sid in genera[genus]:
            local = []
            for offset, r in enumerate(species[sid]):
                if r['path'] not in cache:
                    batch = [v for v in species[sid][offset:offset + 16] if v['path'] not in cache]
                    with ThreadPoolExecutor(max_workers=32) as pool:
                        for v, q in zip(batch, pool.map(inspect, [v['path'] for v in batch])): cache[v['path']] = q
                q = cache[r['path']]; why = None
                if q['fingerprint'] in fingerprints: why = 'previous_sequence'
                elif q['fingerprint'] in selected_fp or q['fingerprint'] in {v['fingerprint'] for v in group + local}: why = 'duplicate_sequence'
                elif r['biosample'] in selected_samples or r['biosample'] in {v['biosample'] for v in group + local}: why = 'duplicate_biosample'
                elif q['acgt_fraction'] < .99 or q['n50'] < 10000: why = 'sequence_qc'
                if why: failed.append(dict(object_id=r['object_id'], reason=why)); continue
                local.append(dict(r, **q))
                if len(local) == CONFIG['genomes_per_species']: break
            if len(local) == CONFIG['genomes_per_species']:
                group.extend(local); group_species.append(sid)
            if len(group_species) == CONFIG['species_per_genus']: break
        atomic(cachepath, cache)
        if len(group_species) != CONFIG['species_per_genus']: continue
        chosen.append(genus); selected.extend(group)
        selected_samples.update(r['biosample'] for r in group)
        selected_fp.update(r['fingerprint'] for r in group)
        phy[group[0]['phylum']] += 1; fam[group[0]['family']] += 1
        print('Selected unseen genera', len(chosen), 'genomes', len(selected), flush=True)
    if len(chosen) < 10: raise RuntimeError('Too few independent genera; cannot support planned validation')
    selected.sort(key=lambda r: r['object_id'])
    assert not ({accession_root(r['object_id']) for r in selected} & used)
    assert not (selected_samples & samples or selected_fp & fingerprints or set(chosen) & old_genera)
    assert len(selected_fp) == len(selected) == len(selected_samples)
    tsv(DATA / 'validation_genomes.tsv', selected)
    if failed: tsv(DATA / 'qc_exclusions.tsv', failed)
    metadata = [BASE / 'as.bacteria.txt'] + [BASE / (r + '_real_bacteria_speciesTaxid_origanismName') for r in lineage]
    atomic(complete, dict(genomes=len(selected), species=len(selected) // 5, genera=len(chosen),
        requested_genera=CONFIG['validation_genera'], families=len(fam), phyla=len(phy),
        prior_accession_roots=len(used), prior_biosamples=len(samples), prior_fingerprints=len(fingerprints),
        overlap_accessions=0, overlap_biosamples=0, overlap_fingerprints=0, overlap_genera=0,
        qc_rejections=len(failed), metadata_rejections=dict(rejected),
        prior_manifest_sha256=manifests, catalog_sha256={str(p): sha(p) for p in metadata},
        protocol_sha256=sha(HERE / 'PROTOCOL.md'), config_sha256=sha(HERE / 'config.json')))
    return selected


def sample_segments(records, length, number, rng):
    slots = np.array([max(0, len(s) - length + 1) for s in records], dtype=np.int64)
    if not slots.sum(): raise ValueError('no contig long enough for requested read')
    cumulative = np.cumsum(slots)
    flat = rng.integers(0, int(cumulative[-1]), number)
    indices = np.searchsorted(cumulative, flat, side='right')
    offsets = flat - np.r_[0, cumulative[:-1]][indices]
    return [records[int(i)][int(p):int(p) + length] for i, p in zip(indices, offsets)]


def derive(records, donor, condition, source_id):
    method = condition['method']
    # Error-rate variants share positions and base-level error uniforms.
    seed_identity = (source_id, method, condition.get('length'), condition.get('coverage'), condition.get('fraction'))
    seed = int(key(seed_identity)[:16], 16); rng = np.random.default_rng(seed)
    provenance = dict(seed=seed, method=method, requested=condition)
    if method == 'fragment':
        length = min(condition['length'], max(map(len, records)))
        result = sample_segments(records, length, 1, rng)
        provenance['truncated'] = length != condition['length']
    elif method == 'tiles':
        length = condition['tile_length']
        tiles = [s[i:i + length] for s in records for i in range(0, len(s) - length + 1, length)]
        count = min(len(tiles), max(1, round(sum(map(len, records)) * condition['fraction'] / length)))
        result = [tiles[int(i)] for i in rng.choice(len(tiles), count, replace=False)]
        contam = round(len(result) * condition['contamination'])
        if contam:
            donor_tiles = [s[i:i + length] for s in donor for i in range(0, len(s) - length + 1, length)]
            if contam > len(donor_tiles): raise ValueError('insufficient unique donor tiles')
            result[:contam] = [donor_tiles[int(i)] for i in rng.choice(len(donor_tiles), contam, replace=False)]
        provenance.update(contaminant_records=contam, total_records=len(result))
    elif method == 'reads':
        length = condition['length']
        count = max(1, round(sum(map(len, records)) * condition['coverage'] / length))
        result = sample_segments(records, length, count, rng)
        rate = condition['substitution']; mutations = 0
        if rate:
            a = np.frombuffer(b''.join(result), dtype=np.uint8).copy()
            letters = np.frombuffer(b'ACGT', dtype=np.uint8)
            lookup = np.zeros(256, dtype=np.uint8); lookup[letters] = np.arange(4)
            valid = np.isin(a, letters)
            mask = (rng.random(len(a)) < rate) & valid
            changes = rng.integers(1, 4, len(a), dtype=np.uint8)
            a[mask] = letters[(lookup[a[mask]] + changes[mask]) % 4]
            raw = a.tobytes(); result = [raw[i:i + length] for i in range(0, len(raw), length)]
            mutations = int(mask.sum())
        provenance.update(substitutions=mutations, total_records=len(result))
    else: raise ValueError('unexpected derived method')
    return result, provenance


def materialize_source(task):
    source, donor, split = task
    records = fasta_records(source['path']); donor_records = fasta_records(donor['path'])
    result = []
    for condition in CONFIG['conditions']:
        if condition['method'] == 'full': continue
        oid = source['object_id'] + '__' + condition['name']
        path = DATA / 'queries' / split / (oid + '.fa')
        meta = path.with_suffix('.json')
        if meta.exists():
            row = json.loads(meta.read_text())
            if row['sha256'] != sha(path): raise ValueError('derived query changed')
        else:
            seqs, provenance = derive(records, donor_records, condition, source['object_id'])
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('wb') as f:
                for i, s in enumerate(seqs): f.write(b'>' + str(i).encode() + b'\n' + s + b'\n')
            lengths = list(map(len, seqs))
            row = dict(object_id=oid, path=str(path), source=source['object_id'], donor=donor['object_id'],
                split=split, genus=source['genus'], species=source['species'], kind=condition['kind'],
                condition=condition['name'], bases=sum(lengths), records=len(seqs),
                median_record_length=float(np.median(lengths)), sha256=sha(path),
                source_sha256=sha(source['path']), donor_sha256=sha(donor['path']), provenance=provenance)
            atomic(meta, row)
        result.append(row)
    return result


def prepare_split(split, genomes):
    completion = DATA / (split + '_preparation.json')
    if completion.exists(): return
    grouped = collections.defaultdict(list)
    for r in genomes: grouped[(r['genus'], r['species'])].append(r)
    source_rows, reference_rows, tasks = [], [], []
    for _, group in sorted(grouped.items()):
        group = sorted(group, key=lambda r: key(('role', r['object_id'])))
        reference_rows.extend(group[:CONFIG['references_per_species']])
        source_rows.extend(group[CONFIG['references_per_species']:])
    for source in source_rows:
        donors = [r for r in reference_rows if r['genus'] == source['genus'] and r['species'] != source['species']]
        donor = min(donors, key=lambda r: key(('donor', source['object_id'], r['object_id'])))
        tasks.append((source, donor, split))
    objects = []
    refids = {r['object_id'] for r in reference_rows}
    for r in genomes:
        records = fasta_records(r['path']); lengths = list(map(len, records))
        objects.append(dict(object_id=r['object_id'], path=r['path'], source=r['object_id'], donor='',
            split=split, genus=r['genus'], species=r['species'], kind='assembly', condition='full',
            bases=sum(lengths), records=len(lengths), median_record_length=float(np.median(lengths)),
            sha256=sha(r['path']), role='reference' if r['object_id'] in refids else 'query'))
    print('Deriving', split, 'queries from', len(tasks), 'source genomes', flush=True)
    with ThreadPoolExecutor(max_workers=16) as pool:
        for i, rows in enumerate(pool.map(materialize_source, tasks), 1):
            objects.extend(dict(r, role='query') for r in rows)
            if i % 30 == 0: print(split, 'derived source genomes', i, '/', len(tasks), flush=True)
    metadata = [{k: v for k, v in r.items() if k != 'provenance'} for r in objects]
    tsv(DATA / (split + '_objects.tsv'), metadata)
    pairs, workloads = {}, []
    def add_pair(a, b):
        if a == b: raise ValueError('self pair')
        pid = key((a, b))[:24]
        if (a, b) not in pairs: pairs[a, b] = dict(pair_id=pid, a=a, b=b)
        return pid
    for genus in sorted({r['genus'] for r in genomes}, key=key):
        refs = sorted(r['object_id'] for r in reference_rows if r['genus'] == genus)
        full = sorted(r['object_id'] for r in genomes if r['genus'] == genus)
        for i, a in enumerate(full):
            for b in full[i + 1:]: add_pair(a, b)
        for condition in CONFIG['conditions']:
            queries = sorted(r['object_id'] for r in objects if r['genus'] == genus and
                             r.get('role') == 'query' and r['condition'] == condition['name'])
            if len(queries) != 9 or len(refs) != 6: raise ValueError('unbalanced query/reference cohort')
            retrieval = [add_pair(q, r) for q in queries for r in refs]
            cluster = [add_pair(a, b) for i, a in enumerate(queries) for b in queries[i + 1:]]
            workloads.append(dict(workload_id=f'{split}_{genus}_{condition["name"]}', split=split,
                genus=genus, condition=condition['name'], kind=condition['kind'],
                queries=queries, references=refs, retrieval_pairs=retrieval, cluster_pairs=cluster))
    tsv(DATA / (split + '_pairs.tsv'), list(pairs.values()))
    atomic(DATA / (split + '_workloads.json'), workloads)
    atomic(completion, dict(genomes=len(genomes), genera=len(genomes) // 15, objects=len(objects),
        pairs=len(pairs), workloads=len(workloads), derived_query_files=len(objects) - len(genomes),
        query_sources=len(source_rows), reference_genomes=len(reference_rows),
        simulated_reads=True, numpy_version=np.__version__, config_sha256=sha(HERE / 'config.json'),
        files={str(DATA / (split + s)): sha(DATA / (split + s))
               for s in ['_objects.tsv', '_pairs.tsv', '_workloads.json']}))
    print('Prepared', split, json.loads(completion.read_text()), flush=True)


def main():
    DATA.mkdir(parents=True, exist_ok=True); RUN.mkdir(parents=True, exist_ok=True)
    identity = dict(protocol_sha256=sha(HERE / 'PROTOCOL.md'), config_sha256=sha(HERE / 'config.json'),
                    prepare_sha256=sha(__file__), common_sha256=sha(HERE / 'scenario_common.py'))
    path = DATA / 'registration.json'
    if path.exists() and json.loads(path.read_text()) != identity:
        raise ValueError('registered preparation code changed; document amendment before rerun')
    atomic(path, identity)
    validation = select_validation()
    development = read_tsv(OLD / 'natural_candidates.tsv')
    tsv(DATA / 'development_genomes.tsv', development)
    prepare_split('development', development)
    prepare_split('validation', validation)
    metadata = {}
    for split in ['development', 'validation']:
        for r in read_tsv(DATA / (split + '_objects.tsv')): metadata[r['path']] = r['sha256']
    atomic(RUN / 'input_sha256.json', metadata)
    atomic(DATA / 'complete.json', dict(registration=identity, input_manifest_sha256=sha(RUN / 'input_sha256.json')))


if __name__ == '__main__': main()
