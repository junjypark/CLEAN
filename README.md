# CLEAN

**Clusterwise inference leveraging spatial autocorrelation in neuroimaging**

R code to apply CLEAN to neuroimaging data. The current version supports parallel computing using the *doParallel* package.

* Park, J.Y., Fiecas, M. CLEAN: Leveraging spatial autocorrelation in neuroimaging data in clusterwise inference.

## Background
CLEAN currently supports group-level inference for neuroimaging data registered in the cortical surface. Please refer [R: ciftiTools](https://github.com/mandymejia/ciftiTools), [Python: ciftify](https://github.com/edickie/ciftify), or [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) for processing cortical surface data.

## Installation
To install the latest development builds directly from GitHub, please run the followings:

```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("junjypark/CLEAN")
```

Then the package can be loaded directly in R:
```R
library(CLEAN)
```

## Usage
- [One-sample test](#id-section1)
- [Two-sample test](#id-section2)
- [General linear model (GLM)](#id-section3)
- [Miscellaneous](#id-section4)

<div id='id-section1'/>

### One-sample test (test for the grand mean)

<div id='id-section2'/>

### Two-sample test (test for the group difference in means)

<div id='id-section3'/>

### General linear model (GLM)

<div id='id-section4'/>

### Miscellaneous
Please refer the [SpLoc package](https://github.com/junjypark/SpLoc), which conducts clusterwise inference for longitudinal neuroimaging data in comparing two groups's growth/decay. It currently does not support leveraging spatial autocorrelations.


## Questions?
Please forward your inquiries to **junjy.park [[at]] utoronto [[dot]] ca**.
