
# **B16 In Vivo Screen**

This repository supplies the code developed in the study of (insert authorship) **_"B16 In Vivo Screen"_**. The corresponding pre-print can be found on bioRxiv ( doi: ). It can be used to reproduce the results of the study and investigate the methodology to be used for other datasets.

## **Requirements**

* R version 3.6.1.
* R libraries: stats, colorspace, MASS, EBImage, flowViz, coda

## **Data**

A set of example data from the paper is provided in the script B16 In Vivo Screen Example MCMC Analysis. This data represents a group of C57BL/6 mice injected with B16F0 wild type cells that was used as a control group for comparison with the DNMT3A and PTPN11 knock outs injected in C57BL/6 mice.

All raw data used in the experiment is provided in the B16_InVivoScreen_RawData.R script.

## **Quick start**

To reproduce the results, download the script. Running the scripts will generate all relevant figures and data tables for the given portion of the study.

# General notes

The code provided in this repository reproduces the main results of the study of (insert authors) **_"B16 In Vivo Screen"_** but it is not meant as a self-contained module for use in further analysis.

For a single example of the MCMC analysis, see the B16 In Vivo Screen Example MCMC Analysis script.

For the full analysis of all data used in the paper, first run the B16_InVivoScreen_RawData.R script, then run the B16_InVivoScreen_MCMC.R script. Note that the pdf generating lines are currently commented out to speed up the code run time. To see the generated figures, simply uncomment these lines.

To generate the summary figure used in the paper, complete the full analysis of all data, then run the ComparisonFigure.R file.

## Citation

