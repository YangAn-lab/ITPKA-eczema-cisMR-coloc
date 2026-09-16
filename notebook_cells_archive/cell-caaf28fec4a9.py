# Confirm Ensembl ID mapping for the unnamed S2 entries
for eid in ['ENSG00000243708', 'ENSG00000260814', 'ENSG00000247556', 'ENSG00000103966']:
    rr = requests.get(f'https://rest.ensembl.org/lookup/id/{eid}?content-type=application/json', timeout=30)
    j = rr.json()
    print(eid, '->', j.get('external_name'), j.get('biotype'), f"chr{j.get('seq_region_name')}:{j.get('start')}-{j.get('end')}")

# check PLA2G4B / JMJD7 ids in our candidate table
print(cand[cand['external_name'].isin(['PLA2G4B', 'JMJD7', 'NUSAP1', 'OIP5', 'EHD4'])][['id', 'external_name']].to_string(index=False))