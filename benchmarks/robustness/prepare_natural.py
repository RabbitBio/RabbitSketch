#!/usr/bin/env python3
"""Taxonomy/QC selection and strict exclusion, before any sketch outcome."""
from pathlib import Path
import collections,csv,gzip,hashlib,json,re,time
from concurrent.futures import ThreadPoolExecutor
from common import *
BASE=Path('/home/bigssd/zt/ref')

def inspect(path):
    raw=Path(path).read_bytes()
    records=fasta_records(gzip.decompress(raw) if str(path).endswith('.gz') else raw)
    lengths=sorted((len(s) for _,s in records),reverse=True);total=sum(lengths);running=0
    for length in lengths:
        running+=length
        if running>=total/2:n50=length;break
    valid=sum(len(s)-s.count(b'N') for _,s in records)
    return dict(fingerprint=fingerprint(records),source_sha256=hashlib.sha256(raw).hexdigest(),bases=total,n50=n50,acgt_fraction=valid/total,records=len(records))

def main():
    if (DATA/'natural_preparation.json').exists():print('Natural preparation already complete');return
    metadata=BASE/'as.bacteria.txt'
    with metadata.open() as f:
        for line in f:
            if line.startswith('#assembly_accession'):fields=line.lstrip('#').rstrip().split('\t');break
        rows=list(csv.DictReader(f,fieldnames=fields,delimiter='\t'))
    prior_paths={};used=set();prior_manifests={}
    for manifest in sorted((ROOT/'runs').glob('*/input_sha256.json')):
        if manifest.parent.name=='robustness_20261001':continue
        prior_manifests[str(manifest)]=sha(manifest)
        for path in json.loads(manifest.read_text()):
            root=accroot(path)
            if root:
                used.add(root)
                if root not in prior_paths and Path(path).exists():prior_paths[root]=path
    for r in rows:
        a,b=accroot(r['assembly_accession']),accroot(r['gbrs_paired_asm'])
        if a in used or b in used:
            used.add(a)
            if b:used.add(b)
    used_samples={r['biosample'] for r in rows if accroot(r['assembly_accession']) in used and r['biosample'] not in ['','na']}
    frozen=DATA/'prior_fingerprints.json'
    if frozen.exists():prior_info=json.loads(frozen.read_text())
    else:
        ordered=sorted(prior_paths.values());print('Fingerprinting earlier genomes',len(ordered),flush=True)
        prior_info={}
        with ThreadPoolExecutor(max_workers=64) as pool:
            for i,(path,result) in enumerate(zip(ordered,pool.map(inspect,ordered))):
                prior_info[path]=result
                if (i+1)%1000==0:print('prior genomes',i+1,flush=True)
        atomic(frozen,prior_info)
    used_fingerprints={v['fingerprint'] for v in prior_info.values()}
    files={p.name.split('_',2)[0]+'_'+p.name.split('_',2)[1]:p for p in (BASE/'bacteria').glob('GC[AF]*.fna.gz')}
    lineage={}
    for rank in ['genus','family','phylum']:
        path=BASE/(rank+'_real_bacteria_speciesTaxid_origanismName')
        lineage[rank]={accroot(r['assembly_accession']):r['genus_id'] for r in read_tsv(path)}
    species=collections.defaultdict(list);rejected=collections.Counter()
    for r in rows:
        root=accroot(r['assembly_accession']);acc=r['assembly_accession']
        if acc not in files:rejected['missing_local_version']+=1;continue
        if root in used:rejected['earlier_accession']+=1;continue
        if r['biosample'] in used_samples:rejected['earlier_biosample']+=1;continue
        if r['biosample'] in ['','na']:rejected['missing_biosample']+=1;continue
        if r['group']!='bacteria' or r['genome_rep']!='Full' or r['excluded_from_refseq'] not in ['','na']:rejected['metadata_quality']+=1;continue
        if not (1_000_000<=int(r['genome_size'])<=10_000_000 and int(r['contig_count'])<=300):rejected['metadata_size_contigs']+=1;continue
        if any(root not in lineage[x] or lineage[x][root] in ['0','1','2','na',''] for x in lineage):rejected['missing_lineage']+=1;continue
        rr=dict(object_id=acc,path=str(files[acc]),species=r['species_taxid'],biosample=r['biosample'],organism=r['organism_name'],assembly_level=r['assembly_level'],genus=lineage['genus'][root],family=lineage['family'][root],phylum=lineage['phylum'][root])
        species[rr['species']].append(rr)
    genera=collections.defaultdict(list)
    for sid,rs in species.items():
        rs.sort(key=lambda r:key(r['object_id']));bybio={}
        for r in rs:bybio.setdefault(r['biosample'],r)
        if len(bybio)>=5:
            species[sid]=list(bybio.values());genera[rs[0]['genus']].append(sid)
    genera={g:sorted(ss,key=key) for g,ss in genera.items() if len(ss)>=3}
    print('Eligible genera with >=3 species:',len(genera),flush=True)
    cachepath=DATA/'candidate_qc.json';cache=json.loads(cachepath.read_text()) if cachepath.exists() else {}
    selected=[];selected_samples=set();selected_fp=set();phylum_counts=collections.Counter();family_counts=collections.Counter();qc_fail=[]
    remaining=set(genera);chosen_genera=[]
    while remaining and len(chosen_genera)<60:
        def order(g):
            r=species[genera[g][0]][0];return phylum_counts[r['phylum']],family_counts[r['family']],key(g)
        genus=min(remaining,key=order);remaining.remove(genus);genus_objects=[];good_species=[]
        for sid in genera[genus]:
            opts=species[sid]
            chosen=[];local_fp=set();local_bio=set()
            for offset,r in enumerate(opts):
                # Read the next deterministic batch only when needed. This is the
                # same first-five-passing selection without reading thousands of
                # unused assemblies from an overrepresented species.
                if r['path'] not in cache:
                    new=[x for x in opts[offset:offset+16] if x['path'] not in cache]
                    with ThreadPoolExecutor(max_workers=64) as pool:
                        for x,v in zip(new,pool.map(inspect,[x['path'] for x in new])):cache[x['path']]=v
                q=cache[r['path']];reason=None
                if q['fingerprint'] in used_fingerprints:reason='earlier_sequence'
                elif q['fingerprint'] in selected_fp|local_fp or q['fingerprint'] in {x['fingerprint'] for x in genus_objects}:reason='duplicate_sequence'
                elif r['biosample'] in selected_samples|local_bio or r['biosample'] in {x['biosample'] for x in genus_objects}:reason='duplicate_biosample'
                elif q['acgt_fraction']<.99 or q['n50']<10000:reason='sequence_qc'
                if reason:qc_fail.append(dict(object_id=r['object_id'],reason=reason));continue
                chosen.append(dict(r,**q));local_fp.add(q['fingerprint']);local_bio.add(r['biosample'])
                if len(chosen)==5:break
            if len(chosen)==5:genus_objects+=chosen;good_species.append(sid)
            if len(good_species)==3:break
        atomic(cachepath,cache)
        if len(good_species)<3:continue
        selected+=genus_objects;selected_samples.update(r['biosample'] for r in genus_objects);selected_fp.update(r['fingerprint'] for r in genus_objects);chosen_genera.append(genus)
        phylum_counts[genus_objects[0]['phylum']]+=1;family_counts[genus_objects[0]['family']]+=1
        print(f'genera {len(chosen_genera)}/60; genomes {len(selected)}',flush=True)
    if len(chosen_genera)<40:raise RuntimeError('fewer than 120 QC-passing species')
    outgen=DATA/'natural_genomes';outgen.mkdir(exist_ok=True)
    def materialize(r):
        dst=outgen/(r['object_id']+'.fna')
        if not dst.exists():dst.write_bytes(gzip.decompress(Path(r['path']).read_bytes()))
        return dict(r,source_path=r['path'],path=str(dst),sha256=sha(dst))
    with ThreadPoolExecutor(max_workers=64) as pool:selected=list(pool.map(materialize,selected))
    selected.sort(key=lambda r:r['object_id']);tsv(DATA/'natural_candidates.tsv',selected)
    tsv(DATA/'natural_candidate_objects.tsv',[dict(object_id=r['object_id'],path=r['path']) for r in selected])
    bygen=collections.defaultdict(list)
    for r in selected:bygen[r['genus']].append(r)
    pairs=[]
    def pair(a,b,kind):
        a,b=sorted([a,b],key=lambda r:r['object_id'])
        return dict(pair_id=a['object_id']+'__'+b['object_id'],a=a['object_id'],b=b['object_id'],pair_class=kind,species_a=a['species'],species_b=b['species'],genus_a=a['genus'],genus_b=b['genus'],family_a=a['family'],family_b=b['family'])
    for rs in bygen.values():
        for i,a in enumerate(rs):
            for b in rs[i+1:]:pairs.append(pair(a,b,'within_species' if a['species']==b['species'] else 'within_genus'))
    representatives=[min(rs,key=lambda r:key(r['object_id'])) for rs in bygen.values()]
    between={'cross_genus_same_family':[],'cross_family':[]}
    for i,a in enumerate(representatives):
        for b in representatives[i+1:]:
            kind='cross_genus_same_family' if a['family']==b['family'] else 'cross_family'
            between[kind].append(pair(a,b,kind))
    # Up to 120 controls of each class, with deficit reported rather than mislabeled.
    control_counts={}
    for kind,rs in between.items():
        rs.sort(key=lambda r:key(r['pair_id']));pairs+=rs[:120];control_counts[kind]=min(120,len(rs))
    tsv(DATA/'natural_candidate_pairs.tsv',pairs)
    if qc_fail:tsv(DATA/'natural_qc_exclusions.tsv',qc_fail)
    (DATA/'natural_fastani.txt').write_text(''.join(r['path']+'\n' for r in selected))
    tax_ids={r[k] for r in selected for k in ['species','genus','family','phylum']};names={}
    with (BASE/'names.dmp').open() as f:
        for line in f:
            fields=[x.strip() for x in line.split('|')]
            if fields[0] in tax_ids and fields[3]=='scientific name':names[fields[0]]=fields[1]
    atomic(DATA/'taxonomy_names.json',names)
    info=dict(genomes=len(selected),species=len({r['species'] for r in selected}),genera=len(chosen_genera),families=len(family_counts),phyla=len(phylum_counts),candidate_pairs=len(pairs),controls=control_counts,prior_accession_roots=len(used),prior_biosamples=len(used_samples),prior_sequence_fingerprints=len(used_fingerprints),metadata_exclusions=dict(rejected),qc_exclusions=len(qc_fail),source_metadata={str(p):sha(p) for p in [metadata,BASE/'names.dmp',*[BASE/(r+'_real_bacteria_speciesTaxid_origanismName') for r in lineage]]},prior_manifest_hashes=prior_manifests,protocol_sha256=sha(HERE/'PROTOCOL.md'))
    assert not ({accroot(r['object_id']) for r in selected}&used)
    assert not ({r['biosample'] for r in selected}&used_samples)
    assert not ({r['fingerprint'] for r in selected}&used_fingerprints)
    assert len({r['fingerprint'] for r in selected})==len(selected)
    atomic(DATA/'natural_preparation.json',info);print(json.dumps(info,indent=2),flush=True)

if __name__=='__main__':main()
