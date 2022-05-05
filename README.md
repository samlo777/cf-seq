# cf-seq

### Introduction: 

CF-Seq is an easy to use web-application geared toward microbiologists who study cystic fibrosis pathogens. The app allows its users to see the broad set of treatment conditions that their pathogen(s) of interest have been exposed to in prior experiments, filter studies by these conditions, and view  analysis results (volcano plots, MA plots, tables of p values and fold changes) for differentially expressed genes. Users are able to zero in on expression of individual genes, and also see how biological pathways are up or down regulated in certain conditions  

### Important Notes: 

If you use CF-Seq for your own research, please cite the following publication: 

CF-Seq, An Accessible Web Application for Rapid Re-Analysis of Cystic Fibrosis Pathogen RNA Sequencing Studies
Samuel L. Neff, Thomas H. Hampton, Charles Puerner, Liviu Cengher, Georgia Doing, Alexandra J. Lee, Katja Koeppen, Ambrose L. Cheung, Deborah A. Hogan, Robert A. Cramer, Bruce A. Stanton
bioRxiv 2022.03.07.483313; doi: https://doi.org/10.1101/2022.03.07.483313

If you have any questions about the application or its source code, please reach out to neff [dot] sam1 [at] gmail [dot] com

### Requirements: 

You are welcome to simply access the application online, no downloads from GitHub required: http://scangeo.dartmouth.edu/CFSeq/

But given that you are here in Github, you probably wish to download the app to your own computer. It requires 30 GB of disk space to run. You also need to have R (https://www.r-project.org/) and the RStudio IDE (https://www.rstudio.com/products/rstudio/download/) installed.

When you download the directory from Github (hit the green ‘Code’ button up above and then ‘Download ZIP’), you will need to do a couple things before the app will run. First, open up ‘CFSeqV1.RProj’. Once opened in RStudio, you can select the files ‘app.R’ and ‘Data Setup.R’. Next, you need to execute the code as specified in ‘Data Setup.R’ - this will load in the study data on which the app is based. Finally, with the app.R script opened, you can hit ‘Run App’ to see the app in action. It may take a minute to boot up because loading the data into the app requires some time.
