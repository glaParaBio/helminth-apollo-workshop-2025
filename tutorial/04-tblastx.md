# Preparing a synteny track

In this part you will make a synteny track for your own mitochondrial genomes
in order to assess the quality of the assembly and annotation.

The example code here uses `tblastx` to find regions of sequence similarity
(i.e., putative synteny) between the mitochondrial genomes of the two
*Trichuris* species used so far. Replace these with your own genomes. Note that
`tblastx` searches translated nucleotide database against translated nucleotide
query, so it is quite sensitive but slow on larger genomes. A faster but less
sensitive alternative could use `minimap2`.

First prepare a database for one of the two species:

```
makeblastdb -in data/trichuris_trichiura/trichuris_trichiura.MITO.fa -dbtype nucl
```

Now perform the blast search. We opt for tabular output with columns suitable
for later steps:

```sh
mkdir -p tblastx

HEADER='qaccver qlen qstart qend sstrand saccver slen sstart send nident'
echo "#${HEADER}" | sed 's/ /\t/g' > tblastx/TMUE_MITO_vs_TTRE_MITO.out

tblastx -query data/trichuris_muris/trichuris_muris.MITO.fa \
    -db data/trichuris_trichiura/trichuris_trichiura.MITO.fa \
    -evalue 0.1 \
    -max_target_seqs 10 \
    -outfmt "6 $HEADER" >> tblastx/TMUE_MITO_vs_TTRE_MITO.out
```

At the moment JBrowse2/Apollo3 cannot read tblastx output directly so we will
convert it to [paf](https://lh3.github.io/minimap2/minimap2.html#10) format
with the script [blast2paf.py](../scripts/blast2paf.py) provided in this
repository. Note that we write the output to the `jbrowse_data` directory since
this is the directory visible to JBrowse:

```sh
./scripts/blast2paf.py > jbrowse_data/TMUE_MITO_vs_TTRE_MITO.paf
```

## Loading to Apollo

At this point you should be able to load the synteny track in Apollo together
with your genomes and annotations using similar commands as in
[02-loading-data](02-loading-data.md). 

However, since mitochondrial genomes have an atypical genetic code, the checks
built in Apollo may give false positives about missing stop codons, internal
stop codons, etc. So, after loading your assembly, but before importing
features, issues this command (replace the assembly name with you own):

```
apollo assembly check --delete --check CDSCheck --assembly trichuris_muris.MITO
```

Using the example data files these are the commands to load assemblies and
synteny (edit as appropriate for your files):

```
apollo assembly add-from-fasta -e data/trichuris_muris/trichuris_muris.MITO.fa -a trichuris_muris.MITO --force 
apollo assembly add-from-fasta -e data/trichuris_trichiura/trichuris_trichiura.MITO.fa -a trichuris_trichiura.MITO --force

apollo assembly check -d -c CDSCheck -a trichuris_muris.MITO
apollo feature import -d data/trichuris_muris/trichuris_muris.MITO.gff3 -a trichuris_muris.MITO

MURIS_MITO_ID=$(
  apollo assembly get |
    jq --raw-output '.[] | select(.name=="trichuris_muris.MITO")._id'
)
TRICHIURA_MITO_ID=$(
  apollo assembly get |
    jq --raw-output '.[] | select(.name=="trichuris_trichiura.MITO")._id'
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
