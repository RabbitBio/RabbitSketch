#!/usr/bin/env python3
"""Fetch and audit small published benchmark assets, without downloading genomes."""
import argparse
import csv
import datetime as dt
import json
from pathlib import Path
import re
import shutil
import tarfile
import urllib.request

from run import DEFAULT_ROOT, atomic_json, file_sha


def fetch(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        tmp = path.with_suffix(path.suffix + '.part')
        with urllib.request.urlopen(url, timeout=30) as src, tmp.open('wb') as dest:
            shutil.copyfileobj(src, dest)
        tmp.replace(path)
    return dict(url=url, path=str(path), bytes=path.stat().st_size, sha256=file_sha(path))


def mash(root):
    root.mkdir(parents=True, exist_ok=True)
    archive = root / 'Escherichia.tar.gz'
    source = fetch('https://gembox.cbcb.umd.edu/mash/Escherichia.tar.gz', archive)
    members = ['README.txt', 'ids.txt', 's500.txt', 's1000.txt', 's5000.txt', 's10000.txt']
    assets = []
    # Only copy allowlisted ordinary files; never extract archive paths/links.
    with tarfile.open(archive) as tar:
        for name in members:
            member = tar.getmember('Escherichia/' + name)
            if not member.isfile():
                raise ValueError('unexpected archive member type')
            dest = root / name
            with tar.extractfile(member) as src, dest.open('wb') as out:
                shutil.copyfileobj(src, out)
            assets.append(dict(path=name, bytes=dest.stat().st_size, sha256=file_sha(dest)))
    genomes = []
    for line in (root / 'ids.txt').read_text().splitlines():
        identifier, original_name = line.split('\t', 1)
        assembly = re.search(r'GC[AF]_\d+\.\d+', original_name)
        genomes.append(dict(genome_id=identifier, original_name=original_name,
                            assembly_accession=assembly.group(0) if assembly else '',
                            retrieval_status='not_downloaded' if assembly else 'needs_accession_resolution'))
    if len(genomes) != 500 or len({g['genome_id'] for g in genomes}) != 500:
        raise ValueError('Mash ID map does not contain 500 unique genomes')
    with (root / 'genome_manifest.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(genomes[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(genomes)
    valid_ids = {g['genome_id'] for g in genomes}
    tables, baseline = {}, None
    for size in (500, 1000, 5000, 10000):
        pairs, observations = {}, []
        conflicting_pairs, conflicting_ani = set(), set()
        examples = []
        with (root / f's{size}.txt').open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            columns = reader.fieldnames
            for line_number, row in enumerate(reader, 2):
                a, b = row['id'].split('_')
                key = tuple(sorted((int(a), int(b))))
                if a not in valid_ids or b not in valid_ids or a == b:
                    raise ValueError('invalid Mash pair')
                ani = float(row['nuc'])
                if not 0 <= ani <= 1 or any(not 0 <= float(row[c]) <= 1 for c in columns[2:]):
                    raise ValueError('similarity outside [0,1]')
                observations.append((line_number, a, b, ani))
                if key in pairs and row != pairs[key]:
                    conflicting_pairs.add(key)
                    if row['nuc'] != pairs[key]['nuc']:
                        conflicting_ani.add(key)
                    if len(examples) < 3:
                        examples.append(dict(pair_id=row['id'], first=pairs[key], later=row,
                                             later_line_number=line_number))
                else:
                    pairs[key] = row
        if baseline is None:
            baseline = observations
        elif observations != baseline:
            raise ValueError('ANI observations/pair IDs differ between Mash tables')
        tables[str(size)] = dict(unique_pairs=len(pairs), rows=len(observations),
                                duplicate_rows=len(observations)-len(pairs),
                                conflicting_pairs=len(conflicting_pairs),
                                conflicting_ani_pairs=len(conflicting_ani), conflict_examples=examples,
                                k_values=[int(c[1:]) for c in columns[2:]],
                                jaccard_truth_type='finite_sketch_estimate_not_exact')
    # Preserve both observations with their original row identity, never choose
    # an arbitrary duplicate or average conflicting values into a fake truth.
    with (root / 'ani_observations_unresolved.tsv').open('w') as stream:
        stream.write('source_line\tgenome_a\tgenome_b\tani_observation\n')
        for line, a, b, ani in baseline:
            stream.write(f'{line}\t{a}\t{b}\t{ani:.6f}\n')
    audit = dict(source=source, extracted_assets=assets, genomes=len(genomes),
                 assembly_accessions=sum(bool(g['assembly_accession']) for g in genomes),
                 pairs_expected=500*499//2, pairs_observed=tables['1000']['unique_pairs'], tables=tables,
                 fasta_included=False, ani_method='Nucmer/dnadiff (archive README)',
                 exact_jaccard_available=False,
                 ani_reference_usable=all(not t['duplicate_rows'] for t in tables.values()),
                 status='requires_pair_id_reconciliation' if any(t['duplicate_rows'] for t in tables.values()) else 'audited',
                 next_step='Resolve and retrieve original assembly versions, checksum FASTA, compute exact canonical k-mer sets for each k. Do not substitute current assembly versions silently.')
    atomic_json(root / 'audit.json', audit)
    return {k: audit[k] for k in ('genomes', 'assembly_accessions', 'pairs_observed', 'fasta_included', 'exact_jaccard_available', 'status')}


def dashing(root):
    manifest_path = root / 'snapshot.json'
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        commit = previous['commit']
    else:
        url = 'https://api.github.com/repos/langmead-lab/dashing-experiments/commits/master'
        request = urllib.request.Request(url, headers={'User-Agent': 'RabbitSketch2-benchmark'})
        with urllib.request.urlopen(request, timeout=30) as response:
            commit = json.load(response)['sha']
    names = ['README.md', 'LICENSE', 'accuracy/pairselector.py', 'accuracy/pairwise_benchmark.cpp',
             'accuracy/genomes_for_exp.txt', 'dsexp/dsexp.cpp']
    assets = []
    for name in names:
        url = f'https://raw.githubusercontent.com/langmead-lab/dashing-experiments/{commit}/{name}'
        assets.append(fetch(url, root / name))
    genomes, pairs = {}, []
    for number, line in enumerate((root / 'accuracy/genomes_for_exp.txt').read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        paths = line.split()
        if len(paths) != 2:
            raise ValueError('Dashing pair manifest must have two paths per line')
        accessions = []
        for path in paths:
            found = re.search(r'GC[AF]_\d+\.\d+', Path(path).name)
            if not found:
                raise ValueError('Dashing path has no versioned assembly accession')
            accession = found.group(0)
            if accession in genomes and genomes[accession] != path:
                raise ValueError('conflicting paths for one Dashing accession')
            genomes[accession] = path
            accessions.append(accession)
        pairs.append(dict(source_line=number, assembly_a=accessions[0], assembly_b=accessions[1]))
    with (root / 'pairs.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['source_line', 'assembly_a', 'assembly_b'], delimiter='\t')
        writer.writeheader(); writer.writerows(pairs)
    with (root / 'genomes.tsv').open('w') as stream:
        stream.write('assembly_accession\toriginal_path\tretrieval_status\n')
        for accession, path in sorted(genomes.items()):
            stream.write(f'{accession}\t{path}\tnot_downloaded\n')
    atomic_json(manifest_path, dict(commit=commit, assets=assets,
                fetched_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                pairs=len(pairs), unique_genomes=len(genomes), fasta_downloaded=False,
                role='published_methodology_and_pair_manifest; not RabbitSketch results',
                license='See downloaded LICENSE; reference files stay outside the RabbitSketch source tree'))
    return dict(commit=commit, files=len(assets), pairs=len(pairs), unique_genomes=len(genomes))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=DEFAULT_ROOT / 'sources')
    p.add_argument('--dataset', choices=['mash', 'dashing', 'all'], default='all')
    args = p.parse_args()
    if args.dataset in ('mash', 'all'):
        print('Mash:', json.dumps(mash(args.root / 'mash500')), flush=True)
    if args.dataset in ('dashing', 'all'):
        print('Dashing:', json.dumps(dashing(args.root / 'dashing')), flush=True)


if __name__ == '__main__':
    main()
