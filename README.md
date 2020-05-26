# articSlurmy

## Introduction
articSlurmy is a [Bash](https://www.gnu.org/software/bash/) based implementation of the [ARTIC networks nCoV-2019 novel coronavirus bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) for [Slurm](https://slurm.schedmd.com/documentation.html) clusters, with added downstream analysis and QC. It is not dependent on a workflow system such as Nextflow and submits Slurm job arrays for each [ONT](https://nanoporetech.com/) flowcell of output directly.  articSlurmy was developed at [Northumbria University Newcastle](https://www.northumbria.ac.uk/), and is used by [NU-OMICS](https://www.northumbria.ac.uk/business-services/engage-with-us/research/nu-omics/) on their HPC facility for their [COG-UK](https://www.cogconsortium.uk/) sequencing effort. A typical run of a [GridION](https://nanoporetech.com/products/gridion) should only take about 10-30 minutes to analyse depending on depth of data produced.


## Features

* Submit script with input/output sanity checking
* Upload ready output for [COG-UK](https://www.cogconsortium.uk/)
* Extensive QC with, global depth, N counts and amplicon level readouts, using [mosdpeth](https://github.com/brentp/mosdepth) and [bedtools](https://bedtools.readthedocs.io/en/latest/)
* Variant annotation with [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
* Lineage assignment with [Pangolin](https://github.com/hCoV-2019/pangolin)


## Setting up the enviroments
Clone this repo into your home dir:

`git clone https://github.com/MattBashton/articSlurmy`

Next setup the artic network conda environment as [explained here](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html):

```
git clone https://github.com/artic-network/artic-ncov2019.git
cd artic-ncov2019
conda env create -f environment.yml
cd
```

Next activate this environment and install additional channels and dependencies:

```
conda activate artic-ncov2019
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install pigz
conda install mosdepth
conda install bedtools
```

Next we need to install [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/) for variant annotation, which sadly is not available for easy conda based installation, via the [registration form](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/download/), articSlurmy will assume you have installed it to `~/annovar` with the "SARS-CoV-2" annotation package provided by the author (see "2020Apr28" update on main page) into `~/annovar/sarscov2db`

Finally [pangolin](https://github.com/hCoV-2019/pangolin) (Phylogenetic Assignment of Named Global Outbreak LINeages) need to be installed into its own conda environment to provide lineage assignment:

```
conda deactivate
cd
git clone https://github.com/hCoV-2019/pangolin
cd pangolin
conda env create -f environment.yml
conda activate pangolin
pip install .
conda deactivate
```

Finally you should edit/check some local config settings in `artic.sh`, specifically check the `# Some defaults` section which details a few paths and run settings related to amplicons and primers you might want to check. I find `$TMPDIR` is often not configured correctly on Slurm clusters so create my own temp dir using `mktemp`, however the `-p` prefix will need pointing at the right path on the worker nodes. _e_._g_. `/scratch`, `/tmp`, or `/local` _etc_. depending on your local setup.

## Usage

Before submitting your first run you will need to create a tab delimited `.tsv` file, format: `barcode  central_sample_id` - this ensures output is in ready to upload format and created with correct `central_sample_id` prefix. This should look like:

```
barcode01	CENT-188D9E
barcode02	CENT-188F98
barcode03	CENT-18994E
```

Submit a run via:

`submit_ARTIC_run.sh <tab delimited barcode dir - sample ID mapping> <input run dir> <sequencing_summary.txt> <output dir for analysis> <global run name>`

_e_._g_.:

`submit_ARTIC_run.sh sample_sheet.tsv  run_output/ run_output/sequencing_summary.txt analysis_output 2020-05-25_MACHINE-ID_CENT-0001_FLOWCELL-ID`

## Ouput

The output directory specified above will include two subdirs, i) `processed/` which contains all output generated for all samples in the run, prefixed with their `central_sample_ids`, and ii) `upload/` which within a global run name directory has subdirectories for each `central_sample_id` that contain `alignment.bam` and `consensus.fa` ready for submission to CLIMB if you are contributing to COG-UK:

```
upload/
└── 2020-05-25_MACHINE-ID_CENT-0001_FLOWCELL-ID/
    ├── CENT-188D9E
    │   ├── alignment.bam
    │   └── consensus.fa
    ├── CENT-188F98
    │   ├── alignment.bam
    │   └── consensus.fa
    ├── CENT-18994E
    │   ├── alignment.bam
    │   └── consensus.fa
```

The output directory will also contain `artic.oJOBID.N` and `artic.eJOBID.N` files in Sun Grid Engine style, the standard out files contained detailed QC data as outlined in Features above, in addition to files written to the `processed/` directory.

## Run stats report
Once a run is complete the `getStats.sh` script located in the output dir can be run to generate a run report which is also written to the file `run_stats.tsv` this contains details of total reads, depth, consensus Ns, variants called and lineage assignment per barcode/sample.  This takes the form of:

```
central_sample_id  barcode    total_reads  post_guppyplex_reads  reads_carried_into_analysis  input_depth  ouput_depth  number_consensus_Ns  %N     variants_called  lineage
CENT-188D9E        barcode01  540824       532110                98.39                        6723.14      401.12       124                  0.41   14               B.1.1
CENT-188F98        barcode02  179183       174903                97.61                        2200.65      371.30       658                  2.20   3                B.3
CENT-18994E        barcode03  3026802      2837018               93.73                        35952.53     431.82       123                  0.41   11               B.1.56
```

## Similar work
There are two Nextflow based pipelines which more or less have similar functionality to this pipeline, the [Connar lab ncov2019-artic-nf](https://github.com/connor-lab/ncov2019-artic-nf) and the [nf-core artic](https://github.com/nf-core/artic) implementation.
