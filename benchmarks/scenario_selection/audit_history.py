#!/usr/bin/env python3
"""Distinguish held-out genera from taxa never seen anywhere in project history."""
from scenario_common import *


def main():
    used,manifests=set(),{}
    for pattern in ['*/input_sha256.json','*/real_input_sha256.json']:
        for p in sorted((ROOT/'runs').glob(pattern)):
            if p.parent==RUN:continue
            manifests[str(p)]=sha(p)
            used.update(filter(None,(accession_root(x) for x in json.loads(p.read_text()))))
    path=BASE/'genus_real_bacteria_speciesTaxid_origanismName'
    lineage={accession_root(r['assembly_accession']):r['genus_id'] for r in read_tsv(path)}
    prior={lineage[a] for a in used if a in lineage and lineage[a] not in ['0','1','2','na','']}
    validation=read_tsv(DATA/'validation_genomes.tsv');development=read_tsv(DATA/'development_genomes.tsv')
    vg={r['genus'] for r in validation};dg={r['genus'] for r in development}
    result=dict(prior_direct_accession_roots=len(used),prior_roots_with_genus_mapping=sum(a in lineage for a in used),
        genera_in_all_prior_input_manifests=len(prior),development_genera=len(dg),validation_genera=len(vg),
        validation_genus_overlap_with_development=len(vg&dg),validation_genus_overlap_with_all_prior_inputs=len(vg&prior),
        validation_genera_unseen_in_all_prior_inputs=len(vg-prior),overlap_ids=sorted(vg&prior),
        lineage_sha256=sha(path),prior_manifest_sha256=manifests,code_sha256=sha(__file__),
        note='Genus novelty is relative to the 60-genus development panel. Accession, BioSample and sequence exclusions apply to all historical inputs.')
    atomic(OUT/'historical_genus_audit.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['prior_manifest_sha256','overlap_ids']},indent=2))


if __name__=='__main__':main()
