#!/usr/bin/env python3
"""Trace data exclusions, taxonomy coverage and original source snapshots."""
import collections,datetime
from common import *

def main():
    prep=json.loads((DATA/'natural_preparation.json').read_text());nat=read_tsv(DATA/'natural_candidates.tsv')
    used=set()
    for manifest in prep['prior_manifest_hashes']:
        for path in json.loads(Path(manifest).read_text()):
            root=accroot(path)
            if root:used.add(root)
    accuracy_paths=read_tsv(ROOT/'data/regimes_20261001/ani_objects.tsv')
    accuracy_roots={accroot(r['path']) for r in accuracy_paths if accroot(r['path'])}
    metadata=Path('/home/bigssd/zt/ref/as.bacteria.txt')
    with metadata.open() as f:
        for line in f:
            if line.startswith('#assembly_accession'):fields=line.lstrip('#').rstrip().split('\t');break
        old_species=set();old_samples=set();accuracy_species=set()
        for r in csv.DictReader(f,fieldnames=fields,delimiter='\t'):
            if accroot(r['assembly_accession']) in used or accroot(r['gbrs_paired_asm']) in used:
                old_species.add(r['species_taxid']);old_samples.add(r['biosample'])
            if accroot(r['assembly_accession']) in accuracy_roots or accroot(r['gbrs_paired_asm']) in accuracy_roots:accuracy_species.add(r['species_taxid'])
    old_genera=set();accuracy_genera=set();lineage=Path('/home/bigssd/zt/ref/genus_real_bacteria_speciesTaxid_origanismName')
    for r in read_tsv(lineage):
        if accroot(r['assembly_accession']) in used:old_genera.add(r['genus_id'])
        if accroot(r['assembly_accession']) in accuracy_roots:accuracy_genera.add(r['genus_id'])
    prior_fp={v['fingerprint'] for v in json.loads((DATA/'prior_fingerprints.json').read_text()).values()}
    overlaps=dict(accession_roots=len({accroot(r['object_id']) for r in nat}&used),biosamples=len({r['biosample'] for r in nat}&old_samples),sequence_fingerprints=len({r['fingerprint'] for r in nat}&prior_fp),species=len({r['species'] for r in nat}&old_species),genera=len({r['genus'] for r in nat}&old_genera))
    assert not any(overlaps[k] for k in ['accession_roots','biosamples','sequence_fingerprints'])
    names=json.loads((DATA/'taxonomy_names.json').read_text());species={r['species'] for r in nat};unnamed=[s for s in species if any(w in names.get(s,'').lower() for w in [' sp.','bacterium','uncultured','unclassified'])]
    sources={path:dict(sha256=value,filesystem_mtime_utc=datetime.datetime.fromtimestamp(Path(path).stat().st_mtime,datetime.timezone.utc).isoformat()) for path,value in prep['source_metadata'].items()}
    # These are snapshot/file dates, not invented download or release dates.
    accuracy_overlap=dict(accession_roots=len({accroot(r['object_id']) for r in nat}&accuracy_roots),species=len(species&accuracy_species),genera=len({r['genus'] for r in nat}&accuracy_genera))
    info=dict(overlap_with_earlier_inputs=overlaps,overlap_with_earlier_accuracy_panel=accuracy_overlap,new_genomes=len(nat),new_species=len(species),new_genera=len({r['genus'] for r in nat}),unnamed_or_placeholder_species=len(unnamed),taxonomic_identities_not_required_to_be_unseen=True,source_snapshots=sources,limitations=['exact sequence fingerprint does not detect circular-origin rotations, alternate contig boundaries or near-identical isolates','BioSample/accession exclusion is limited to available assembly metadata','taxonomic balance is a curated validation design, not an estimate of the RefSeq population distribution','all genomes are bacterial and 1–10 Mb; this does not establish large-eukaryote or read-workflow generalization'])
    info.update(prior_accuracy_species_taxids=sorted(accuracy_species),prior_accuracy_genus_taxids=sorted(accuracy_genera))
    atomic(ANALYSIS/'data_provenance.json',info)
    print(json.dumps(info,indent=2),flush=True)

if __name__=='__main__':main()
