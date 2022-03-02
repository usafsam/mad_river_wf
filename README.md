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
