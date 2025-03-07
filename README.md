# FT-SEM
The simulation and emprical study code of the paper “A robust and powerful GWAS method for family trios supporting within-family Mendelian randomization analysis”


# Installation
It is easy to install the development version of FTSEM package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute. This R package is currently running smoothly on version 4.4.0.
```
# install.packages("devtools")
library(devtools)
install_github("ShunZhang0816/FT-SEM")
```
# Usage
There are two main functions in FTSEM package, one is FT_SEM for performing analysis based on family trios data with one offspring, another one is process_family_data to transform the data from PLINK 1.9 to what the FT_SEM function needs. You can find the instructions by '?FT_SEM' and '?process_family_data'.
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
The first five columns of this file represent individual information: FID (Family ID), IID (Individual ID), PAT (Paternal ID), MAT (Maternal ID), and y (Phenotypic value). The subsequent columns correspond to SNP genotypes, where each SNP is coded as 0, 1, or 2, indicating the number of minor alleles.  
Another straightforward example demonstrating the use of FT-SEM to conduct GWAS separately on exposure and outcome datasets, followed by combining the summary results using IVW to estimate causal effects, can be found in ./Example/Example_2.  
Example_2 consists of two files, Example_2_Exposure and Example_2_Outcome, both of which follow the same format as Example_1.

## Expected output
Example_1 and Example_2 are executed on a system powered by an AMD 7500F CPU and an RX 7900 GRE GPU, running the Windows 11 operating system. Upon running Example_1, you will receive a data.frame that includes point estimates, confidence intervals, and p-values. Similarly, running Example_2 will also generate a data.frame containing point estimates, confidence intervals, and p-values


# Results reproduced
All results for all methods used in the FT-SEM paper can be reproduced at ./Simulation and ./Empirical_Study. It is important to note that reproducing the empirical analysis requires obtaining publicly available data in advance. For details, please refer to the links provided in the article.
