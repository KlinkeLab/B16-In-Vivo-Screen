
# **B16 In Vivo Screen**

This repository supplies the code developed in the study of Pirkey et al. **_"Head-to-head comparison of CCN4, DNMT3A, PTPN11, and SPARC as suppressors of anti-tumor immunity"_**. The corresponding pre-print can be found on bioRxiv ( doi:https://doi.org/10.1101/2022.04.01.486749 ) and has now been published as Pirkey et al. (2023) Cell Mol Bioeng Oct 28; 16(5-6):431-442. It can be used to reproduce the results of the study and investigate the methodology to be used for other datasets.

## **Requirements**

* R version 3.6.1.
* R libraries: stats, colorspace, MASS, EBImage, flowViz, coda

## **Data**

All raw data used in the experiment is provided in the B16_InVivoScreen_RawData.R script.

## **Quick start**

To reproduce the results, download the scripts. Running the scripts in the following order will generate all relevant figures and data tables for the given portion of the study.

1) 1_B16_InVivoScreen_RawData.R
2) 2_B16_InVivoScreen_MCMC.R
3) 3_B16_InVivoScreen_Analysis_Normalized.R

# General notes

The code provided in this repository reproduces the main results of the study of Pirkey et al. **_"Head-to-head comparison of CCN4, DNMT3A, PTPN11, and SPARC as suppressors of anti-tumor immunity"_** but it is not meant as a self-contained module for use in further analysis.

For the full analysis of all data used in the paper, first run the B16_InVivoScreen_RawData.R script, then run the B16_InVivoScreen_MCMC.R script. Note that the pdf generating lines are currently commented out to speed up the code run time. To see the generated figures, simply uncomment these lines.

To generate the summary figure used in the paper, complete the full analysis of all data, then run the ComparisonFigure.R file.

## Citation
Pirkey, A.C.; Deng, W.; Norman, D.; Razazan, A.; Klinke, D.J.; **_"Head-to-head comparison of CCN4, DNMT3A, PTPN11, and SPARC as suppressors of anti-tumor immunity”_**, (2023) Cell Mol Bioeng Oct 28; 16(5-6):431-442.
