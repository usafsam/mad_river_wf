# Mad River workflow

Named after the Mad River, one of the major rivers that runs through Dayton, Ohio.

Mad River is a workflow developed by Padraic Fanning at United States Air Force School of Aerospace Medicine (USAFSAM) Applied Technologies & Genomics Division, based on the [Cecret workflow by Erin Young](https://github.com/UPHL-BioNGS/Cecret) and work by Dr. Anthony Fries.
This workflow is initially designed for SARS-COV-2 sequencing with the Illumina Nextera XT library prep workflow using [version 1 of the "midnight" primer set by Freed and Silander](https://www.protocols.io/view/sars-cov2-genome-sequencing-protocol-1200bp-amplic-btsrnnd6).
Currently, this workflow (in its current state) has been tested on data generated from NextSeq runs and includes potentially helpful diagnostics such as spike gene coverage and variant coverage/quality metrics.
The tools used are mostly sourced from [the Docker images provided by StaPH-B](https://github.com/StaPH-B/docker-builds), which include:

- [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/)
- [iVar](https://andersen-lab.github.io/ivar/html/)
- [nextclade](https://github.com/nextstrain/nextclade) (_no Docker image provided by StaPH-B as of 2021-08-17_)
- [pangolin](https://cov-lineages.org/resources/pangolin.html)
- [samtools](https://www.htslib.org/)
- [VADR](https://github.com/ncbi/vadr)

However, the Docker images and the exact Conda environments can be customized to your liking.

# Getting Started

In order to run this workflow, you will need [Conda](https://docs.conda.io/en/latest/miniconda.html), [Nextflow](https://www.nextflow.io), and [Docker](https://www.docker.com/).
Then, if you do not wish to modify this workflow's scripts or reference files, run:

```none
nextflow run usafsam/mad_river_wf \
    --reads {READS_DIR} \
    --run_info {PATH_TO}/RunInfo.xml \
    --stats_json {PATH_TO}/Stats.json \
    --outdir {OUTDIR}
```

If you are using this workflow locally, replace `usafsam/mad_river_wf` with the path to where this workflow resides.

# Specifying Primers

The default set of primers used in this workflow is version 2.0 of the SARS-CoV-2 Midnight Amplicon panel [as provided by IDT (released April 2022)](https://www.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels/sars-cov-2-midnight-amp-panel).
A TSV file of these primers can be found in the `reference/` directory of this repository, along with prior versions of this set.
A python script, `process_primer_reference.py`, takes this TSV file and produces a FASTA file of the primers (this gets used by BBDuk), along with BED files for the primers and amplicons.
The version numbers follow the principles of [semantic versioning](https://semver.org/), where the MAJOR version corresponds to a new release from IDT and the MINOR version corresponds to one or more additional primers being spiked into the reaction mixture.
Here is a table that shows which primers are included in each version of the set.

| Primer Name | v1.0 | v1.1 | v1.2 | v2.0 |
|---|---|---|---|---|
| 1\_LEFT through 29\_RIGHT | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| 28\_LEFT\_OMICRON (C27807T) | :x: | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| 22\_RIGHT\_OMICRON (G22599A) | :x: | :x: | :white_check_mark: | :white_check_mark: |
| 23\_LEFT\_OMICRON (C22522T) | :x: | :x: | :white_check_mark: | :white_check_mark: |
| 26\_LEFT\_OMICRON (C25708T) | :x: | :x: | :white_check_mark: | :white_check_mark: |
| 21\_RIGHT\_OMICRON (71\_RIGHT from ARTIC v4.1) | :x: | :x: | :x: | :white_check_mark: |

You can adapt the general format of the TSVs found in `reference/` to the set of primers you have.
To use a different suite of primers other than the default, override the value of `params.primer_tsv` in a user-provided config file, and/or specify `--primer_tsv $PRIMER_SCHEME` when invoking Nextflow from the command line.
