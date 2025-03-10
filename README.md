# FT-SEM
The simulation and empirical study code for the paper entitled “A robust and powerful GWAS method for family trios supporting within-family Mendelian randomization analysis”


# Installation
It is easy to install the development version of FTSEM package using the 'devtools' package. The typical installation time on a "normal" desktop computer is less than one minute. This R package is currently running smoothly on version 4.4.0.
```
# install.packages("devtools")
library(devtools)
install_github("ShunZhang0816/FT-SEM")
```
# Usage
There are two main functions in FTSEM package, one is FT_SEM for performing analysis based on family trio data with one offspring, another one is process_family_data to transform the data from PLINK 1.9 to what the FT_SEM function needs. You can find the instructions by '?FT_SEM' and '?process_family_data'.
```
library(FTSEM)

?FT_SEM

?process_family_data
```

# Example
One simple example to use the FT-SEM to perform GWAS can be found at /Example/Example_1.  
Example_1 is in the standard format extracted using PLINK 1.9. By applying the --recodeA option, traditional PLINK file formats can be converted into this format.  
For example,  
```
plink --bfile your_data_name --recodeA --out output_data_name
```
After that, you will get a file with a .raw extension. You can directly import it into R for use or change its extension to .txt before importing it into R.  
The first five columns of this file represent individual information: FID (Family ID), IID (Individual ID), PAT (Paternal ID), MAT (Maternal ID), and PHENOTYPE (Phenotypic value). The subsequent columns correspond to SNP genotypes, where each SNP is coded as 0, 1, or 2, indicating the number of minor alleles.  
  
Another straightforward example demonstrating the use of FT-SEM to conduct GWAS separately for exposure and outcome datasets, followed by combining the summary results using IVW to estimate causal effects, can be found in ./Example/Example_2.  
Example_2 consists of two data files, Example_2_Exposure and Example_2_Outcome, both of which have the same format as Example_1.

## Expected output
Example_1 and Example_2 are executed on a system powered by an AMD 7500F CPU and an RX 7900 GRE GPU, running the Windows 11 operating system. Upon running Example_1, you will receive a data.frame that includes point estimates, confidence intervals, and p-values of association effect sizes. Similarly, running Example_2 will also generate a data.frame containing point estimates, confidence intervals, and p-values of causal effect sizes.

Example_1 output

| SNP  | Method | Beta_o_e           | SE_o_e             | p_wald_o_e           | CI_lower_o_e       | CI_upper_o_e      | Beta_f_e           | SE_f_e             | p_wald_f_e        | CI_lower_f_e        | CI_upper_f_e      | Beta_m_e            | SE_m_e             | p_wald_m_e        | CI_lower_m_e       | CI_upper_m_e      |
| ---- | ------ | ------------------ | ------------------ | -------------------- | ------------------ | ----------------- | ------------------ | ------------------ | ----------------- | ------------------- | ----------------- | ------------------- | ------------------ | ----------------- | ------------------ | ----------------- |
| SNP1 | FT-SEM | 0.175176278737423  | 0.0497620176361808 | 0.000431090037002146 | 0.0776427241705088 | 0.272709833304337 | 0.0247095927763713 | 0.0553931836387655 | 0.655542048660935 | -0.0838610471556091 | 0.133280232708352 | -0.0108437983896527 | 0.0570975397594303 | 0.849374103157615 | -0.122754976318136 | 0.101067379538831 |
| SNP2 | FT-SEM | 0.0140173158674012 | 0.0643153076412046 | 0.827470557822518    | -0.11204068710936  | 0.140075318844162 | 0.104000057566962  | 0.0713433547406661 | 0.144912238916003 | -0.0358329177247435 | 0.243833032858668 | 0.0735158730036606  | 0.0737678440443646 | 0.318966344324958 | -0.071069101323294 | 0.218100847330615 |
| .... |        |                    |                    |                      |                    |                   |                    |                    |                   |                     |                   |                     |                    |                   |                    |                   |

Example_2 output

| Beta_IVW  | SE_IVW     | CI_Lower   | CI_Upper  | P_Value     |
| --------- | ---------- | ---------- | --------- | ----------- |
| 0.1149939 | 0.03701849 | 0.04243771 | 0.1875502 | 0.001893854 |


# Results reproduced
All results for all methods used in the FT-SEM paper can be reproduced at ./Simulation and ./Empirical_Study. It is important to note that reproducing the empirical analysis requires obtaining publicly available data in advance. For details, please refer to the links provided in in the "Data availability" section of the article.

# Reference
Mc, N. et al. OpenMx 2.0: Extended Structural Equation and Statistical Modeling. Psychometrika 81, 535–549 (2016).
