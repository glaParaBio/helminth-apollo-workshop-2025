
```
makeblastdb -in data/trichuris_trichiura.MITO.fa -dbtype nucl
```

```
mkdir -p tblastx

header='qaccver qlen qstart qend sstrand saccver slen sstart send nident'
echo "#${header}" | sed 's/ /\t/g' > tblastx/TMUE_MITO_vs_TTRE_MITO.out
tblastx -query data/trichuris_muris.MITO.fa \
    -db data/trichuris_trichiura.MITO.fa \
    -evalue 0.1 \
    -max_target_seqs 10 \
    -outfmt "6 $header" >> tblastx/TMUE_MITO_vs_TTRE_MITO.out
cat tblastx/TMUE_MITO_vs_TTRE_MITO.out | scripts/blast2paf.py > jbrowse_data/TMUE_MITO_vs_TTRE_MITO.paf
```

```
apollo assembly add-from-fasta -e data/trichuris_muris.MITO.fa --force
apollo assembly add-from-fasta -e data/trichuris_trichiura.MITO.fa --force

MURIS_MITO_ID=$(
  apollo assembly get |
    jq --raw-output '.[] | select(.name=="trichuris_muris.MITO.fa")._id'
)
TRICHIURA_MITO_ID=$(
  apollo assembly get |
    jq --raw-output '.[] | select(.name=="trichuris_trichiura.MITO.fa")._id'
)

apollo jbrowse get-config > data/config.json
jbrowse add-track \
    data/TMUE_MITO_vs_TTRE_MITO.paf \
  --load inPlace \
  --name "muris_vs_trichiura MITO TBLASTX" \
  --assemblyNames "${MURIS_MITO_ID}","${TRICHIURA_MITO_ID}" \
  --out data/config.json \
  --force

apollo jbrowse set-config data/config.json
rm data/config.json
```
