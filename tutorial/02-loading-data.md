# Loading data

## Setting up the Apollo and JBrowse CLI

Both Apollo and JBrowse have CLI tools that you can use to load data. The
recommended way to install both these tools is through a Node.js package manager
(such as `npm` or `yarn`). See
[the next section](#installing-cli-tools-with-npm-or-yarn) for those
installation instructions.

### Installing CLI tools with `npm` or `yarn`

Use the following to install the Apollo and JBrowse CLI tools.

```bash npm2yarn
npm install -g @apollo-annotation/cli @jbrowse/cli
```

You can test that those commands work by running `apollo --version` and
`jbrowse --version` in your terminal. The output should look something like
this:

```sh-session
$ apollo --version
@apollo-annotation/cli/0.3.5 linux-x64 node-v22.15.0
$ jbrowse --version
@jbrowse/cli/3.4.0 linux-x64 node-v22.15.0
```

If that worked, you can move on to
[configuring the Apollo CLI](#configuring-the-apollo-cli).

## Configuring the Apollo CLI

Open a new terminal in the same directory where you ran the setup commands. To
use the Apollo CLI, we need to configure it with the information for the running
Apollo installation. You can have multiple profiles configured, but we will use
a single default profile. Run these commands:

```sh
apollo config address  http://localhost/apollo
apollo config accessType root
apollo config rootPassword password
apollo login
```

If you need to log in again, run `apollo logout` first, or use
`apollo login --force`.

## Adding assemblies and annotations

### Using the CLI

The next step is to add an assembly. The fasta file of the genome of *Trichuris
trichiura* comes from [WormBase
Parasite](https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS19.genomic_softmasked.fa.gz)
and it has been compressed with `bgzip` and indexed with `samtools faidx`. To add this assembly run:

```sh
apollo assembly add-from-fasta ./data/trichuris_trichiura.PRJEB535.WBPS19.genomic_softmasked.fa.gz \
    --fai ./data/trichuris_trichiura.PRJEB535.WBPS19.genomic_softmasked.fa.gz.fai \
    --gzi ./data/trichuris_trichiura.PRJEB535.WBPS19.genomic_softmasked.fa.gz.gzi \
    --assembly 'Trichuris trichiura' \
    --force
```

Now that we have an assembly, let's add the annotations we want to curate. They
are stored in a GFF3 file. Also the GFF3 file comes from [WormBase
Parasite](https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS19.annotations.gff3.gz).
However to keep the data small we reduced the original file to contain only a
few hundreds of genes. Run this command to import the annotations (it may take
a couple of minutes to complete):

```sh
apollo feature import --delete-existing --assembly 'Trichuris trichiura' ./data/trichuris_trichiura.sample.gff3
```

### Using the graphical interface

Next we're going to add a second assembly and set of annotations. This assembly
is from the related species *Trichuris muris*, also from WormBase. Instead of
the cli, this time we will use the graphical interface. Navigate to
[http://localhost/](http://localhost/), open the *Linear Genome View*. Then
from menu Apollo choose "Add Assembly" (refresh the page if option "Add
Assembly" is not visible). Fill in the submission form with the fasta file and
the two index files for *T. muris* in the `data` directory. Refresh the page once done.

To add features to this assembly use the dialog form "Import Features" from the
Apollo menu. Use the gff file `data/trichuris_muris.sample.gff3`.

The corresponding cli commands would be:

```sh
apollo assembly add-from-fasta ./data/trichuris_muris.PRJEB126.WBPS19.genomic_softmasked.fa.gz \
    --fai ./data/trichuris_muris.PRJEB126.WBPS19.genomic_softmasked.fa.gz.fai \
    --gzi ./data/trichuris_muris.PRJEB126.WBPS19.genomic_softmasked.fa.gz.gzi \
    --assembly 'Trichuris muris' \
    --force

apollo feature import --delete-existing --assembly 'Trichuris muris' ./data/trichuris_muris.sample.gff3
```

## Adding evidence tracks

Apollo is now set up to be able to annotate these genomes. In order to help with
the annotation, though, it's often useful to include evidence tracks. Apollo is
built on top of JBrowse 2, so we'll add these evidence tracks to the underlying
JBrowse configuration. The first thing we need to do is get the IDs of the
assemblies, since we'll need to pass these to JBrowse instead of the Apollo
internal name. You can see the IDs in the output of the assembly adding commands
above, but we'll run the following commands to demonstrate another way to get
them:

```sh
MURIS_ID=$(
  apollo assembly get |
    jq --raw-output '.[] | select(.name=="Trichuris muris")._id'
)
TRICHIURA_ID=$(
  apollo assembly get |
    jq --raw-output '.[] | select(.name=="Trichuris trichiura")._id'
)
```

In order to make the evidence track data available to JBrowse, the
`jbrowse_data/` is visible inside our running application as a directory called
`data/` that is visible to JBrowse. Inside this directory, we have an IsoSeq
file in bam format, as well as a file that shows synteny relationships between
the two assemblies. This particular file was generated with `tblastx`.

The first step is to get the JBrowse configuration stored in Apollo so we can
update it. Run these commands:

```sh
apollo jbrowse get-config > data/config.json
```

Now that we have the configuration, we can use the `jbrowse` CLI tool to add the
evidence tracks. We are using the `inPlace` value for the `--load` flag because
we know how these files are going to be visible in to JBrowse, and we don't want
the CLI to try and copy or alter any files.

```sh
jbrowse add-track \
  data/trichuris_muris_vs_trichuris_trichiura.sample.paf \
  --load inPlace \
  --name "trichuris_muris_vs_trichuris_trichiura TBLASTX" \
  --assemblyNames "${MURIS_ID}","${TRICHIURA_ID}" \
  --out data/config.json \
  --force

jbrowse add-track \
  data/TTRE_all_isoseq.chr2.bam \
  --load inPlace \
  --name "T. trichiura IsoSeq" \
  --assemblyNames "${TRICHIURA_ID}" \
  --out data/config.json \
  --force
```

Now the last step is to send the updated JBrowse config back to Apollo.

```sh
apollo jbrowse set-config data/config.json
rm data/config.json
```

We can now explore these assemblies and start [annotating](03-annotating.md).
