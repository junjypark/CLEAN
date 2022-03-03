# CLEAN

**Clusterwise inference leveraging spatial autocorrelation in neuroimaging**

R code to apply CLEAN to neuroimaging data. The current version supports parallel computing using the *doParallel* package.

> Park, J.Y., Fiecas, M. CLEAN: Leveraging spatial autocorrelation in neuroimaging data in clusterwise inference. [Preprint](https://www.biorxiv.org/content/10.1101/2022.03.02.482664v1)

**Note: This page is currently under construction**


## Contents

1. [Background](#id-background)
2. [Installation](#id-installation)
3. [CLEAN for GLM](#id-cleanglm)
    * Visualization
    * Cautionary notes    
4. [FAQ](#id-tips)
    * [How do I extract surface data from HCP?](#id-q1)    
    * [Which surface should we use for registration?](#id-q2)
    * [How do I obtain a pairwise distance matrix?](#id-q3)
    * [Is it possible to fit CLEAN separately for two hemispheres and combine results afterwards?](#id-q4)
    * [What is the recommended value for max.radius?](#id-q5)
5. [Miscellaneous](#id-misc)
6. [Questions?](#id-question)

---

<div id='id-background'/>

### Background
CLEAN currently supports group-level clusterwise inference for neuroimaging data registered in the cortical surface. Please refer [R: ciftiTools](https://github.com/mandymejia/ciftiTools), [Python: ciftify](https://github.com/edickie/ciftify), or [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) for processing cortical surface data. More helpful information is provided in the [FAQ](#id-tips) section.

<div id='id-installation'/>

---

### Installation
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

---

<div id='id-cleanglm'/>

### CLEAN for GLM (test for the grand mean, test for a difference, general linear model)

Fitting CLEAN consists of three major steps

**Step 1) Obtain new data after leveraging spatial autocorrelation**: This is done by using the spLeverage function with 3 major inputs: (i) a data matrix (ii) a pairwise distance matrix and (iii) covariate information (for two-sample tests or GLM).

For one sample test, use
```R
data.leverage=spLeverage(data, distMat)
```

For two sample test, use
```R
mod0=model.matrix(~1)
data.leverage=spLeverage(data, distMat, mod0)
```

For GLM using potential confounders, use
```R
mod0=model.matrix(~covariates)
data.leverage=spLeverage(data, distMat, mod0)
```
Note: "covariates" above should NOT contain the covariate of interest.


**Step 2) Specify candidate clusters**: Candidate clusters consist of every vertex and its neighbors defined by vertices within a radii. Please use the optional command max.radius from the buildNNmatrixDist_radius() function to specify your neighbors. For example, if you use max.radius=3, then it will create a neighbor information for a vertex, a vertex and its neighbors within 1mm, 2mm, and 3mm.
```R
NNmatrix=buildNNmatrixDist_radius(distMat, max.radius=20)
```

**Step 3) Fit CLEAN**: Once you obtain leveraged data and candidate clusters in Steps 1 and 2, please use Clean() function to obtain the Clean fit.
```R
fit=Clean(data.leverage$out, NNmatrix, seed=NULL)	
```

<div id='id-tips'/>

---

### Frequently asked questions:
<div id='id-q1'/>

**How do I extract surface data from HCP?**
Please refer [ciftiTools](https://github.com/mandymejia/ciftiTools). Once you installed [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) in your computer and obtained data files in nii format and surface information in surf.gii format, then 

```R
library(ciftiTools)
library(rgl)
ciftiTools.setOption("wb_path", "/Applications/workbench")
xii = read_cifti(fname, surfL, surfR, resamp_res = 10242)  #resamp_res: how many vertices to resample
```

Then you can access the cortical data by
```R
xii$data$cortex_left
xii$data$cortex_right
```
and you may collect the data in a matrix format for analysis.The corresponding mesh information can be assessed by
```R
xii$surf$cortex_left
xii$surf$cortex_right
```

<div id='id-q2'/>

**Which surface should we use for registration?**

We use spherical surface as a default and use inflated or midthickness surface for visualization ([reference](https://doi.org/10.1016/j.neuroimage.2016.05.038)). However, please note that there is no definitive answer for this, and there are recent [articles](https://doi.org/10.1016/j.neuroimage.2022.118908) that supported the midthickness surface for registration. Whenever possible, we recommend to conduct an exploratory analysis to make sure the parametric kernel agrees with empirical data. 

<div id='id-q3'/>

**How do I obtain a pairwise distance matrix?**

We recommend using geodesic distance for mesh surfaces. To our knowledge, you may use [Python](https://pypi.org/project/pygeodesic/) or [C++](https://code.google.com/archive/p/geodesic/wikis/ExactGeodesic.wiki) to obtain a pairwise geodesic distance matrix.

<div id='id-q4'/>

**Is it possible to fit CLEAN separately for two hemispheres and combine results afterwards?**: 

Yes, it is necessary to set a brain-wise threshold that controls FWER at the nominal level. Please make sure you specify the same seed and the same number of resamples for the CLEAN() function. Then you may use combine() function to get a new threshold.

```R
Clean.fit.lh=Clean(dataLH, NNmatrixLH, nperm=10000, seed=1)
Clean.fit.rh=Clean(dataRH, NNmatrixRH, nperm=10000, seed=1)
Clean.fit.brain=list(Clean.fit.lh, Clean.fit.rh)
Clean.fit.combine=combine(Clean.fit.brain, alpha=0.05)
```

<div id='id-q5'/>

**What is the recommended value for max.radius?**

The max.radius determines the degree of spatial domain you're borrowing information from. Higher sensitivity obtained from a large value of max.radius, however, comes with the cost of decreased specificity. It should be determined a priori prior to obtaining any result. We empirically found values between 10 and 20 useful for interpretation. 
  
 <div id='id-misc'>

---

### Miscellaneous
If interested, please check out the [SpLoc](https://github.com/junjypark/SpLoc) package, a close family of CLEAN, that conducts clusterwise inference for longitudinal neuroimaging data in comparing two groups's growth/decay. It currently does not support leveraging spatial autocorrelations.

<div id='id-question'>

---

### Questions?
Please forward your inquiries to **junjy.park [[at]] utoronto [[dot]] ca**.

