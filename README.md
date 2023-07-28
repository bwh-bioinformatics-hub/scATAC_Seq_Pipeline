# scATAC_Seq_Pipeline

Functions and Rmarkdowns for QC and analysis of scATAC-seq data.

## Requirements

```
Depends:
  GenomicRanges,
  data.table,
  Matrix
Imports:
  assertthat,
  igraph,
  IRanges,
  Matrix,
  parallel,
  RANN
```

## Installation

This package can be installed from Github using the `remotes` package.

You may first need to register your GitHub PAT, as this is a private repository.
```
Sys.setenv(GITHUB_PAT = "your-access-token-here")
remotes::install_github("bwh-bioinformatics-hub/scATAC_Seq_Pipeline")
```
## Test Data

TBD


To run tests for this package, download the repository from Github and run `remotes::test()` in R from the base directory of the package.

Extra-stringent, CRAN-level package testing can be performed using `remotes::check()` in R.

## Style and Standards

This package aims to conform to the tidyverse style guide:  
https://style.tidyverse.org/index.html

General information about R package conventions can be found in `R Packages`:  
http://r-pkgs.had.co.nz/

### [QC Stage](#stage_qc)

**[Format Arc data: 00_run_arc_formatting.R](#arc_uuid)**
- [Parameters](#arc_uuid_param)
- [Outputs](#arc_uuid_out)

**[Preprocess scATAC Data: 01_crossplatform_preprocessing.sh](#prepro)**
- [Parameters](#prepro_param)
- [Outputs](#prepro_out)

**[Analyze scATAC Quality: 02_run_crossplatform_atac_qc.R](#run_qc)**
- [Parameters](#run_qc_param)
- [Outputs](#run_qc_out)
- [Tests](#run_qc_test)

### [Filter Stage](#stage_filter)

**[Filter fragments based on QC: 03_filter_fragments.sh](#filter)**
- [Parameters](#filter_param)
- [Outputs](#filter_out)
### [Output Stage](#stage_output)

**[Postprocess scATAC Data: 04_crossplatform_postprocessing.sh](#postpro)**
- [Parameters](#postpro_param)
- [Outputs](#postpro_out)

**[Generate output files: 05_run_assemble_atac_outputs.R](#out)**
- [Parameters](#out_param)
- [Outputs](#out_out)
- [Tests](#out_test)

### Pipeline workflow

The scATAC-seq pipeline proceeds through 3 main stages:
- QC 
- Filter 
- Output 

[Return to Contents](#contents)

<a id="depends"></a>

### Setup & Dependencies

This repository requires multiple libraries are installed via `apt-get`:
```
sudo apt-get install -y\
  pandoc \
  libhdf5-devel \
  libxt-dev \
  libcairo2-dev \
  tabix \
  bedtools \
  dos2unix
```


SAMTools and BCFTools are also required for BAM and Tabix processing:
```
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar -xf samtools-1.12.tar.bz2
cd samtools-1.12
./configure --prefix=/usr/bin/
make
make install

cd ../

wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2
tar -xf bcftools-1.12.tar.bz2
cd bcftools-1.12
./configure --prefix=/usr/bin/
make
make install
```

It also depends on several R libraries:

`jsonlite`, `rmarkdown` and `optparse` are available from CRAN, and can be installed in R using:
```
install.packages("jsonlite")
install.packages("rmarkdown")
install.packages("optparse")
```

`scATAC_Seq_Pipeline` and `H5MANIPULATOR` Install with:
```
Sys.setenv(GITHUB_PAT = "[your_PAT_here]")
remotes::install_github("bwh-bioinformatics-hub/scATAC_Seq_Pipeline")
remotes::install_github("bwh-bioinformatics-hub/H5MANIPULATOR")
```

`ArchR` is a key component for doublet detection and .arrow file generation:
```
remotes::install_github("GreenleafLab/ArchR", 
                         ref="master", 
                         repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
ArchR::addArchrGenome("hg38")
```


[Return to Contents](#contents)


## QC

<a id = "arc_uuid"></a>

### Format Arc Outputs

Prior to running processing of ATAC or RNA data from `cellranger-arc` output, run `00_run_arc_formatting.R` to add a common UUID and restructure the metadata files from arc for downstream processing. This ensures that we don't end up with different UUIDs assigned in the RNA and ATAC arms of the pipeline. After this step, `cellranger-arc` outputs are compatible with all downstream steps.

<a id="arc_uuid_param"></a>

### Parameters

There are three parameters for this script:  
- `-t`: path to cellranger-arc outs/
- '-d': Sample ID
- `-o`: path for the HTML output generated by the script

An example run is:
```
Rscript --vanilla ${baseDir}/scripts/00_run_arc_formatting.R \
  -t outs/ \
  -d ${SampleID} \ 
  -o arc_formatting_report.html
```

<a id="arc_uuid_out"></a>

### Outputs

This script outputs .csv files to the outs/ directory for downstream processing:
- arc_singlecell.csv: Arc version of the standard 10x ATAC singlecell.csv output
- atac_summary.csv: Arc version of the standard 10x ATAC summary.csv output
- rna_summary.csv: Arc version of the standard 10x RNA metrics_summary.csv

[Return to Contents](#contents)


<a id="prepro"></a>

## Preprocessing 10x Genomics scATAC Data

This repository utilizes `bedtools` and `GNU parallel` to compute overlaps between fragments and reference datasets. It requires the `fragments.tsv.gz` and `singlecell.csv` files generated by `cellranger-atac count` as inputs, and generates multiple tsv files with overlap counts and sparse matrix data for each reference.

[Return to Contents](#contents)

<a id="prepro_param"></a>

### Parameters

There are 9 parameters for preprocessing:

* `-b` Path to bedtools binary location
* `-g` Genome (hg19 or hg38)  
* `-i` Input fragments.tsv.gz
* `-j` N Parallel Jobs
* `-o` Output Directory
* `-q` Generate QC files only (usually TRUE)
* `-s` Input singlecell.csv
* `-t` Temporary Directory
* `-w` Sample ID

An example run is:
```
bash ${baseDir}/scripts/01_crossplatform_preprocessing.sh \
  -b /shared/apps/bedtools2/bin/bedtools \
  -i SampleID/outs/fragments.tsv.gz \
  -g hg38 \
  -j 16 \
  -o SampleID/preprocessed \
  -q TRUE
  -s SampleID/outs/singlecell.csv \
  -t atac_preprocessing_temp \
  -w SampleID
```

[Return to Contents](#contents)

<a id="prepro_out"></a>

### Output files

For outputs, files are prefixed with the SampleID provided as the `-w` parameter.

The following are generated with `-q TRUE` is used:
- Reference Peaks total counts: [SampleID]_peaks_total_counts.tsv.gz  
- TSS Region total counts: [SampleID]_tss_total_counts.tsv.gz  
- Gene Bodies total counts: [SampleID]_gene_bodies_total_counts.tsv.gz  
- 5kb Window counts: [SampleID]_window_5k_counts.tsv.gz
- UMI count frequencies (saturation curve): [SampleID]_saturation_curve.tsv.gz
- Fragment width counts: [SampleID]_fragment_widths.tsv.gz

Additional outputs are also generated using `-q FALSE`:
- Reference Peaks sparse matrix: [SampleID]_peaks_sparse_matrix.tsv.gz  
- TSS Region sparse matrix: [SampleID]_tss_sparse_matrix.tsv.gz  
- Gene Bodies sparse matrix: [SampleID]_gene_bodies_sparse_matrix.tsv.gz  

[Return to Contents](#contents)

<a id="run_qc"></a>

## Performing 10x scATAC-seq QC

After preprocessing, the resulting files can be analyzed along many quality metrics, perform sequencing saturation projection, remove doublets, generate .h5 for specified references and genomic window sizes, and generate an .arrow file for downstream processing.

<a id="run_qc_param"></a>

### Parameters

There are 10 parameters for this script:  
- `-p or --in_pre`: The path to the preprocessed directory generated by [step 01 (parameter -o), above](#prepro_param).  
- `-t or --in_tenx`: The path to the cellranger-atac outs/  
- `-k or --in_key`: The path to SampleSheet.csv  
- `-s or --in_sample`: SampleID`  
- `-g or --genome`: Either hg38 or hg19. Optional - defaults to hg38. 
* `-q or --qc_only`: Generate QC files only (usually TRUE)
- `-r or --references`: (optional) Either `all` (default) or a comma-separated combination of `gene_bodies,great,grr,peaks,tss`.  
- `-x or --window_sizes`: (optional) either `all` (default) or a comma-separated combination of `window_5k,window_20k,window_100k`.  
- `-d or --out_dir`: A directory to use to output the modified .h5 and JSON metrics  
- `-o or --out_html`: A filename to use to output the HTML summary report file  

An example run for all references and window sizes is:
```
Rscript --vanilla \
  ${baseDir}/scripts/02_run_crossplatform_atac_qc.R \
  -p SampleID/preprocessed \
  -t tenx_out/outs/ \
  -k SampleSheet.csv \
  -w SampleID \
  -g hg38 \
  -q TRUE \
  -d atac_qc/ \
  -o atac_qc/SampleID_qc_report.html
```

An example run for a subset of references and window sizes is:
```
Rscript --vanilla \
  ${baseDir}/scripts/02_run_crossplatform_atac_qc.R \
  -p SampleID/preprocessed \
  -t tenx_out/outs/ \
  -k SampleSheet.csv \
  -w SampleID \
  -g hg38 \
  -q TRUE \
  -r gene_bodies \
  -x window_5k \
  -d atac_qc/ \
  -o atac_qc/SampleID_qc_report.html
```

<a id="run_qc_out"></a>

### Output files

Output files will be prefixed with the PoolID and SampleID supplied in [SampleSheet.csv](#sample_sheet).

Files generated when `q TRUE` is used are:
- All barcode metadata: [SampleID]_all_metadata.csv.gz, e.g. SampleID_all_metadata.csv.gz  
- QC filtered barcode metadata:[SampleID]_filtered_metadata.csv.gz, e.g. SampleID_filtered_metadata.csv.gz  
- Saturation projection: [SampleID]_saturation_projection.csv.gz, e.g. SampleID_saturation_projection.csv.gz  
- Fragment width summary: [SampleID]_fragment_width_summary.csv.gz, e.g. SampleID_fragment_width_summary.csv.gz  
- JSON file: [SampleID]_atac_qc_metrics.json, e.g. SampleID_atac_qc_metrics.json  

Additional outputs are also generated using `-q FALSE`:
- ArchR .arrow file: [SampleID]_archr.arrow, e.g. SampleID.arrow
- Gene TSS .h5 file: [SampleID]_tss.h5, e.g. SampleID_tss.h5  
- Gene Bodies .h5 file: [SampleID]_gene_bodies.h5, e.g. SampleID_gene_bodies.h5  
- 5kb Window .h5 file: [SampleID]_window_5k.h5, e.g. SampleID_window_5k.h5  

[Return to Contents](#contents)

<a id="run_qc_test"></a>

### Tests

Test runs can be performed using datasets provided with the `scATAC_Seq_Pipeline` package by excluding parameters other than `-o`.

```
Rscript --vanilla \
  tenx-atacseq-pipeline/02_run_tenx_atac_qc.R \
  -o test_metadata_report.html
```


## Single Cell ATAC-Seq Report

The `03_scatac_qc_script.R` wrapper script renders `scatac_seq_qc_report_generator.Rmd` to create a QC report (html). 

[Return to Contents](#contents)  

<a id="ngs_report_param"></a>

#### Input Parameters

There are 8 parameters for this script:  

* `-s or --in_sample`:  SampleID
* `-m or --in_method`:  A ";"-delimited string of the data streams being 
processed enclosed in quotes, for example "scatac"
  * `scatac` (single cell ATAC)
* `-i or --in_dir`: The input directory containing all results files to process. Must 
include one subdirectory per modality named in `in_method` above. An example file
tree can be seen in the image [below](#ngs_example_directory),  but the required files 
per subdirectory are as follows:  
  * `scatac`: Contains the following well level files generated by `02_run_crossplatform_atac_qc.R` 
    * `all_metadata.csv.gz`, one file per well. Metadata for all cell barcodes.
    * `filtered_metadata.csv.gz`, one file per well. Metadata for all cell barcodes passing ATAC QC filters.
    * `fragment_width_summary.csv.gz`, one file per well. Fragment width count summary.
    * `saturation_projection.csv.gz`, one file per well. Fragment saturation curve projections.
* `-k or --in_key`: A 6-column .csv Sample Sheet (see format [below](#ngs_sample_sheet)) 
of identifiers for all samples in the batch.
* `-d or --out_dir`: A directory path to use to output the batch report.
* `-o or --out_html`: A filename to use to output the HTML summary report file. 
For example "SampleID_atac_qcreport.html"

Test runs may be performed by supplying only the -b and -o arguments, as seen in
test [below](#ngs_report_test). 

<a id="ngs_example_directory"></a>

[Return to Contents](#contents)  

<a id="ngs_sample_sheet"></a>

##### Sample Sheet Format

`03_scatac_qc_script.R` requires a **Sample Sheet**, provided as the -k parameter.

For example a  run may look like this:
```
TO ADD
```  


Sample sheets should be provided for each dataset. 

An example run:
```
Rscript --vanilla \
    ${baseDir}/scripts/03_scatac_qc_script.R \
    -s SampleID  \
    -m 'scatac' \
    -i ${baseDir}/results/SampleID \
    -k samplesheet   \
    -d ${baseDir}/results/atac_qc \
    -o ${baseDir}/results/atac_qc/SampleID_atac_qcreport.html \

```

[Return to Contents](#contents)

<a id="ngs_report_out"></a>

#### Output Files

`03_scatac_qc_script.R` will generate the HTML reporting file with name as defined by inut parameter -o. In addition the


[Return to Contents](#contents)




## Filtering

<a id="filter"></a>

##Filter fragments.tsv.gz

After QC, we can filter the fragments.tsv.gz file to retain only cell barcodes that pass QC criteria. This reduces the amount of data that goes forward into postprocessing and output.

<a id="filter_param"></a>

### Parameters

There are 6 parameters for this script:  
- `-i`: fragments.tsv.gz file from cellranger-atac outs/  
- `-m`: filtered_metadata.csv.gz file generated in the [previous step](#qc_out)  
- `-j`: N Parallel Jobs  
- `-o`: Path to an output directory  
- `-t`: Path to a temporary directory  
- `-s`: The sample name, e.g. [SampleID]

An example run is:
```
bash ${baseDir}/scripts/04_filter_fragments.sh \
  -i tenx_out/outs/fragments.tsv.gz \
  -m atac_qc/SampleID_filtered_metadata.csv.gz \
  -j 16 \
  -o atac_qc/ \
  -t filtering_temp/ \
  -s SampleID
```

<a id="filter_out"></a>

### Output files

1 filtered fragments.tsg.gz file will be generated. This file will be prefixed with the sample name supplied using the -s parameter:  
- Filtered fragments.tsv.gz: [SampleID]_filtered_fragments.tsv.gz, e.g. SampleID_filtered_fragments.tsv.gz  

This filtered file is compressed using `bgzip` for compatibility with `tabix` indexing.

[Return to Contents](#contents)




<a id="postpro"></a>

### Postprocessing 10x Genomics scATAC Data

This script is similar to preprocessing, but generates a full complement of files for outputs, not just those needed for QC. It utilizes filtered fragments and metadata per well from 04_fragment_filtering and 02_run_crossplatform_atac_qc.

[Return to Contents](#contents)

<a id="postpro_param"></a>

### Parameters

There are 8 parameters for postprocessing:

* `-b` Path to bedtools binary location
* `-g` Genome (hg19 or hg38)  
* `-i` Input filtered_fragments.tsv.gz
* `-j` N Parallel Jobs
* `-m` Input filtered_metadata.csv.gz
* `-o` Output Directory
* `-t` Temporary Directory
* `-s` SampleID

```

An example is:
```
bash ${baseDir}/scripts/05_crossplatform_postprocessing.sh \
  -b /shared/apps/bedtools2/bin/bedtools \
  -i ${baseDir}/results/SampleID/atac_qc/SampleID_filtered_fragments.tsv.gz \
  -g hg38 \
  -j 16 \
  -m ${baseDir}/results/SampleID/atac_qc/SampleID_filtered_metadata.tsv.gz
  -o SampleID/postprocessed \
  -t atac_postprocessing_temp \
  -s SampleID
```

[Return to Contents](#contents)

<a id="postpro_out"></a>

### Output files

For outputs, files are prefixed with the Well or Sample ID provided as the `-w` parameter.

- Reference Peaks total counts: [SampleID]_peaks_total_counts.tsv.gz  
- TSS Region total counts: [SampleID]_tss_total_counts.tsv.gz  
- Gene Bodies total counts: [SampleID]_gene_bodies_total_counts.tsv.gz  
- 5kb Window counts: [SampleID]_window_5k_counts.tsv.gz
- UMI count frequencies (saturation curve): [SampleID]_saturation_curve.tsv.gz
- Fragment width counts: [SampleID]_fragment_widths.tsv.gz
- Reference Peaks sparse matrix: [SampleID]_peaks_sparse_matrix.tsv.gz  
- TSS Region sparse matrix: [SampleID]_tss_sparse_matrix.tsv.gz  
- Gene Bodies sparse matrix: [SampleID]_gene_bodies_sparse_matrix.tsv.gz
## Generating final output files

After postprocessing, the resulting files can be assembled into a set of output files for storage and downstream analysis.

<a id="out_param"></a>

### Parameters

There are 8 parameters (6 required) for this script:  
- `-p or --in_post`: The path to the postprocessed directory generated by [step 05 , above](#postpro_param).  
- `-f or --in_frag`: The path to the filtered_fragments.tsv.gz file  
- `-m or --in_meta`: The path to the filtered_metadata.csv.gz file  
- `-g or --genome`: Either hg38 or hg19. Optional - defaults to hg38. 
- (optional) `-r or --references`:  Either `all` (default) or a comma-separated combination of `gene_bodies,great,grr,peaks,tss`.  
- (optional) `-x or --window_sizes`: Either `all` (default) or a comma-separated combination of `window_5k,window_20k,window_100k`.
- (optional) `-n or --n_cores`: Either `auto` (default) or an integer specifying the number of system cores to use for parallel processing.
- `-d or --out_dir`: A directory to use to output the modified .h5 and JSON metrics  
- `-o or --out_html`: A filename to use to output the HTML summary report file  

An example run foris:
```
Rscript --vanilla \
  tenx-atacseq-pipeline/06_run_assemble_atac_outputs.R \
  -p SampleID/postprocessed \
  -f SampleID/out/SampleID_filtered_fragments.tsv.gz \
  -m SampleID/atac_qc/SampleID_filtered_metadata.csv.gz \
  -g hg38 \
  -d atac_out/ \
  -o atac_out/SampleID_qc_report.html


An example run:
```

Rscript --vanilla \
  tenx-atacseq-pipeline/06_run_assemble_atac_outputs.R \
  -p SampleID/postprocessed \
  -f atac_frag_merged/SampleID_fragments.tsv.gz \
  -m atac_meta_merged/SampleID_metadata.csv.gz \
  -g hg38 \
  -d SampleID/atac_out/ \
  -o SampleID/atac_out/SampleID_out_report.html

```
Output files will be assigned the same prefix used for the postprocessing step, above (auto-detected from filenames).

Files generated when `q TRUE` is used are:
- QC filtered barcode metadata: [SampleID]_filtered_metadata.csv.gz, e.g. SampleID_filtered_metadata.csv.gz  
- Saturation projection: [SampleID]_saturation_projection.csv.gz, e.g. SampleID_saturation_projection.csv.gz  
- Fragment width summary: [SampleID]_fragment_width_summary.csv.gz, e.g. SampleID_fragment_width_summary.csv.gz  
- JSON file: [SampleID]_atac_qc_metrics.json, e.g. SampleID_atac_qc_metrics.json  
- ArchR .arrow file: [SampleID]_archr.arrow, e.g. SampleID.arrow
- Gene TSS .h5 file: [SampleID]_tss.h5, e.g. SampleID_tss.h5  
- Gene Bodies .h5 file: [SampleID]_gene_bodies.h5, e.g. SampleID_gene_bodies.h5  
- 5kb Window .h5 file: [SampleID]_window_5k.h5, e.g. SampleID_window_5k.h5  

