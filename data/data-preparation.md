
```
for genome in trichuris_muris.PRJEB126.WBPS19 trichuris_trichiura.PRJEB535.WBPS19
do
    species=`echo ${genome} | sed 's/\..*//'`
    prj=`echo ${genome} | sed 's/.*PRJ/PRJ/ ; s/\..*//'`
    curl https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/${species}/${prj}/${genome}.genomic_softmasked.fa.gz \
    | gunzip \
    | bgzip > ${genome}.genomic_softmasked.fa.gz
    samtools faidx ${genome}.genomic_softmasked.fa.gz

    ## Also prepare gff3 files. Make them smaller and remove non-CDS records with the same ID
    curl https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/${species}/${prj}/${genome}.annotations.gff3.gz \
    | zcat \
    | grep -P '\tCDS\t|\texon\t|\tgene\t|\tmRNA\t|\tprotein_coding_primary_transcript\t|\trRNA\t|\tsnRNA\t|\ttRNA\t|\tnontranslating_CDS\t|^##' > ${genome}.annotations.gff3
done

awk '$0 ~ "^#" || $1 == "TMUE_LG2" && $4 > 17200000 && $5 < 20001000' trichuris_muris.PRJEB126.WBPS19.annotations.gff3 > trichuris_muris.sample.gff3
awk '$0 ~ "^#" || $1 == "TTRE_chr2" && $4 > 10695000 && $5 < 13024000' trichuris_trichiura.PRJEB535.WBPS19.annotations.gff3 > trichuris_trichiura.sample.gff3

rm trichuris_trichiura.PRJEB535.WBPS19.annotations.gff3 trichuris_muris.PRJEB126.WBPS19.annotations.gff3 
```

```
samtools faidx data/trichuris_trichiura.PRJEB535.WBPS19.genomic_softmasked.fa.gz TTRE_MITO > data/trichuris_trichiura.MITO.fa
samtools faidx data/trichuris_muris.PRJEB126.WBPS19.genomic_softmasked.fa.gz TMUE_MITO > data/trichuris_muris.MITO.fa
```
