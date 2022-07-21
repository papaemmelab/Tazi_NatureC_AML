[![DOI](https://zenodo.org/badge/458019383.svg)](https://zenodo.org/badge/latestdoi/458019383)
# AML_Repo
Github Repository for Tazi et al. Unified classification and risk-stratification in Acute Myeloid Leukemia. 
Please check our online tool : https://www.aml-risk-model.com/
![alt text](data/readme_helper.png?raw=true "Title")

## Installation :

The code was written in Python 3.7 and R 4.1.0.

Make sure to create a dedicated environment as follow :

```
conda env create --name <YOUR_ENV_NAME> 
```

Main packages to download:

```
some_packages <- c("ggplot2","reshape2","data.table","stringr","gridExtra","survival","survminer","bshazard","colorspace","mstate","ggrepel","cmprsk")
for (package in some_packages){
  if(!require(package,character.only = TRUE)){
    install.packages(package,character.only = TRUE)
  }
}
```

## Conclusion :
Using comprehensive molecular profiling data from 3,653 patients we characterize and validate 16 molecular classes describing 100% of AML patients. Each class represents diverse biological AML subgroups, and is associated with distinct clinical presentation, likelihood of response to induction chemotherapy, risk of relapse and death over time. Secondary AML-2, emerges as the second largest class (24%), associates with high-risk disease, poor prognosis irrespective of flow MRD negativity, and derives significant benefit from transplantation. Guided by class membership we derive a 3-tier risk-stratification score that re-stratifies 26% of patients as compared to standard of care. This results in a unified framework for disease classification and risk-stratification in AML that relies on information from cytogenetics and 32 genes. Last, we develop an open-access patient-tailored clinical decision support tool.


