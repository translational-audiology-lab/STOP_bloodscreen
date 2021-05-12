This has the information and description of the code here, which were used to
generate results from given data. This is written in Markdown (`.md`). So, it looks better in a Markdown viewer such as [Typora](https://typora.io/).

## Project info

* Title : Blood biomarkers in chronic tinnitus
* Principal Investigator: Christopher Cederroth (<christopher.cederroth@ki.se>)
* Organization : Karolinska Institutet
* NBIS experts : Mun-Gwan Hong (<mungwan.hong@nbis.se>)

### Aim of the project

*Copied from Redmine*

The overarching aim of this project is to establish novel blood biomarkers for tinnitus, which is the first research priority from the BTA. Biomarkers for tinnitus are currently lacking and these are essential in order to i) understand the pathophysiology of tinnitus; ii) stratify patients; iii) provide read-outs for clinical trials. Similarity in the neuropathophysiology between tinnitus and pain led us to hypothesize that inflammation is involved in tinnitus. 
A preliminary analysis on 548 cases of constant tinnitus and 548 age/sex matched controls showed 10% of the proteins had a strong correlation with smoking, and that 50% of proteins were strongly correlated with age and with self-reported hearing ability (Bonferroni adjusted P < 0.05). We performed a linear regression and adjusted for age, sex, BMI, smoking status, sample collection site, and hearing ability. This allowed us to identify 5 proteins with close to significant adjusted p values, namely FGF-21, MCP-4, CXCL9, GDNF, and MCP-1 (Fig. 1). These proteins were found higher in the tinnitus group and did not associate with stress, anxiety or depression, nor temporomandibular joint disorder, headache or hyperacusis. The goal of this project is to expand the analysis to the neurological panel from O'link on the same samples.

* 

* 

--------------------------------------------------------------------------------

## General info

### Header

Every script has a **header** at the top of it. It has a description about the
script as well as basic info. All **input** and **output** file names are listed
in the header.

### Neighbor folders

Those input and output files are supposed to be stored in the following neighbor
folders. Every script is written to be executed on this current folder accessing
those folders using relative path.

        ../data/raw_internal       # raw input data
               /                   # intermediate derived data from raw data
               /raw_external       # data from external sources, e.g. public database
                 
        ../reports                 # all the reports
        ../results/                # main output folder
                  /figures
                  /tables

### R and Rmd scripts

All R **functions** from installed packages except those listed below, are
called with double colons (`::`) to specify source package clearly.

* Packages in `R-core` (e.g. `base`, `stats`, `utils`, ...),
* Packages in `tidyverse` (e.g. `dplyr`, `ggplot2`, `purrr`, ...)
* `kabelExtra`
* `janitor`
* Packages of the function unable to be accessed using `::`
    - `ggfortify`  # `ggplot2::autoplot` for PCA
    - `lme4`        # `predict` for class `merMod`

All R code followed
[the `tidyverse` **R style guide**](https://style.tidyverse.org/syntax.html).

--------------------------------------------------------------------------------

## Restore environment

Different software environment, e.g. different versions of R packages, from the
one during the development of the code, can yield unexpected error message or
generate dissimilar analysis results. To make the results reproducible, the
software environment around the code are saved. It can be restored ahead of
using the code following the steps below. Here, the process was facilitated by
`conda` management system and `renv` R package. The environment handling
software `conda` is assumed to be pre-installed. Please refer to
[Conda](https://docs.conda.io) for questions about the installation.

1. Open Terminal (Mac / Linux) or Command line (Windows)

2. Navigate to the directory where this README.md file is located using the `cd`
command.

3. Run the commands below. It creates the same Conda environment and activate it
(Note : The `project_a` can be chosen as your preference for your project name).

        conda env create -n project_a -f environment.yml
        conda activate project_a

4. Run `R` and execute this command in `R` to restore the same R environment. 

        renv::restore()

Please refer to 
["introduction to `renv`"](https://rstudio.github.io/renv/articles/renv.html), 
for the R environment handling package, `renv`.

The Conda environment was stored in this file.

* `environment.yml`

The R environment was saved in the folder and file below. They are supposed to
be managed only by the `renv` package under R.

* `renv/` 
* `renv.lock`

--------------------------------------------------------------------------------

## Files

### Master script files

Run these lines below in an R console, which will generate all intermediate data
files for data analysis and the final report. Please make sure the input files
are ready in expected folders. The files were listed below the code lines.  

```R
source("master_data_generation.R")   # data generation
rmarkdown::render("report-5797-Tinnitus.Rmd")    # report writing
```



    ../data/raw_internal/202?-??-??/

The dates in the folder names above indicate when the files were transferred to
NBIS. If multiple versions of one file were transferred, please use the one
delivered on that day. 


### Data generation

#### `master_data_generation.R`

The master file for R data generation. This executes following scripts in the
proper order. 

* `gen-s1-clinc.v01.R` : Read clinical info data
* `gen-s1-olink_proteomic.v01.R` : Read Olink proteomic data
* `gen-s1-olink_proteomic.v02.R` : QC and preprocess - details in `report-QC-Olink_proteomic_data.Rmd`

#### Other data generating scripts


### Heavy analyses

#### `master_heavy_analyses.R`

The master script for computationally heavy analyses.


### Report writing

#### `report-5797-Tinnitus.Rmd`

The master R markdown file to create the final report. Individual chapters were
written in separate R markdown files listed below. This master file runs all of
them in the right order after adding the information about the project and NBIS
support. 

* `report-sample_info.Rmd` : About samples and clinical info
* `report-QC-clinical_info.Rmd` : Clinical info table QC
* `report-QC-Olink_proteomic_data.Rmd` : About QC and preprocessing of Olink
proteomic data


#### Extra analyses

These analyses were not included in the main report, `report.Rmd`. 
Primary reason was that we decided to focus on different analyses.

#### Auxiliary files

* `styles.css` : CSS style, used by Rmarkdown (`.Rmd`) files for reports
* `utils.R` : A collection of useful R functions


### Plots

