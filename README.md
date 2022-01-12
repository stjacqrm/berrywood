# Berrywood
Berrywood is a NextFlow pipeline used to characterize SARS-CoV-2 genomes and provide run-level characterization data.

### Prerequisites

What dependencies you need to run Berrywood


- [NextFlow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)

### Using the pipeline
The pipeline is designed to start from assembled genomes. All fasta files must be in the same directory. Then start the pipeline using:

```
$ nextflow run ~/berrywood/berrywood.nf --assemblies /path/to/assemblies
```

#### To rename the output pdf

```
$ nextflow run ~/berrywood/berrywood.nf --title "Title name in quotes"--assemblies /path/to/assemblies
```

### Author


* **Rachael St. Jacques** - *Bioinformatics Principal Scientist* - [stjacqrm](https://github.com/stjacqrm)
