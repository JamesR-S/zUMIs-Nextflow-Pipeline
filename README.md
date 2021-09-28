# zUMIs Nextflow Pipeline

Ongoing translation of zUMIs Single-Cell RNAseq pipeline into the Nextflow Framework.

Project Also Involves development of a Singularity container containing all dependencies for the pipeline to run.

Part of MSc Bioinformatics dissertation project.

Script can be run with the following command:

```
nextflow run main.nf
```
Nextflow parameters are set in `nextflow.config` file.

A docker container designed to be used for this script can be found at https://hub.docker.com/repository/docker/jamescraufurd/zumis-nf.

Credit to [Parekh et al. 2018](https://doi.org/10.1093/gigascience/giy059) for the original zUMIs pipeline. The pipeline can be found at their [github repository](https://github.com/sdparekh/zUMIs).

