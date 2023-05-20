# CLEAN

**Clusterwise inference leveraging spatial autocorrelation in neuroimaging**

R package to apply CLEAN, CLEAN-R, and CLEAN-V to neuroimaging data. 

### References

**CLEAN**
> Park JY, Fiecas M. (2022) CLEAN: Leveraging spatial autocorrelation in neuroimaging data in clusterwise inference. NeuroImage, 255, 119192. [article link](https://doi.org/10.1016/j.neuroimage.2022.119192)

**CLEAN-R**
> Weinstein SM *et al*. (2022) Spatially-enhanced clusterwise inference for testing and localizing intermodal correspondence. NeuroImage, 264, 119712 [article link](https://doi.org/10.1016/j.neuroimage.2022.119712)

**CLEAN-V**
> Pan R *et al*. (2023+) Spatial-extent inference for testing variance components in reliability and heritability studies. BioRxiv [article link](https://doi.org/10.1101/2023.04.19.537270)


### Update logs

Note: The current version of the package is a beta version and may contain bugs. Please forward your inquiries to `junjy.park [[at]] utoronto [[dot]] ca` and we will respond in 48 hours.

> May 19, 2023: We are currently working on supporting visualizations and preparing manuals. Any big changes will be posted here.

> April 20, 2023: CLEAN-V method implemented in the package.

> Oct 15, 2022: We now support the `get.empirical.variogram' function that computes the group-averaged empirical variogram for exploratory data analyses.

## Contents

1. [Background](#id-background)
2. [Installation](#id-installation)
3. [Usage](#id-usage)
    * [CLEAN: testing activations or group differences via GLM](#id-cleanglm)
    * [CLEAN-R: testing intermodal associations](#id-cleanr)  
    * [CLEAN-V: testing reliability or heritability (forthcoming)](#id-cleanv)
4.  [Visualization](#id-cleanvisualize)
5. [FAQ](#id-tips)
    * [How do I extract surface data from HCP?](#id-q1)    
    * [Which surface should we use for registration?](#id-q2)
    * [How do I obtain a pairwise distance matrix?](#id-q3)
    * [Is it possible to fit CLEAN/CLEAN-R/CLEAN-V separately for two hemispheres and combine results afterwards?](#id-q4)
    * [What is the recommended value for max.radius?](#id-q5)
6. [Miscellaneous](#id-misc)

---

<div id='id-background'/>

### Background
CLEAN supports group-level clusterwise inference for neuroimaging data registered in the cortical surface. Key components of CLEAN include

* an explicit brain-wise spatial covariance modeling of neuroimaging data,
* resampling approaches (sign-flipping or permutation) to control family-wise error rate, and
* a general clusterwise inference for cortical surface.

Compared to classical GLM (or massive-univariate analysis), CLEAN shows superior statistical power. A current implementation is computationally efficient and, using a laptop without parallel computing, takes only a few minutes to analyze 50 subjects' imaging data across 10,000 vertices. The current version also supports parallel computing using the `doParallel` R package.

Please refer [R](https://github.com/mandymejia/ciftiTools), [Python](https://github.com/edickie/ciftify), or [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) for processing cortical surface data. This package is written in R and works well with the [ciftiTools](https://github.com/mandymejia/ciftiTools) R package. More helpful information is provided in the [FAQ](#id-tips) section.

<div id='id-installation'/>

---

### Installation
To install the latest development builds directly from GitHub, please run the followings:

```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("junjypark/CLEAN")
```

Note: If you see an error message, please update your `devtools` package first.

```R
update.packages("devtools")
```

After installation, the package can be loaded directly in R.
```R
library(CLEAN)
```
 

---

<div id='id-cleanglm'/>

### CLEAN: testing activations or group differences via GLM

For one sample test (e.g. testing group-level activation in task-fMRI), use
```R
fit = Clean(ymat = data, distmat = distmat)
```

For two-sample test (e.g. testing difference in means between two groups), use
```R
fit = Clean(ymat = data, group = group, distmat = distmat)
```

For GLM using potential confounders (e.g. brain-behavior associations), use
```R
mod0 = model.matrix(~confounders)
fit = Clean(ymat = data, group = covariate, mod0 = mod0, distmat = distmat)
```
Note: `confounders` above should NOT contain the covariate of interest (`covariate`).

Please refer to the manual for more information about the optional arguments.

```R
help(Clean)
```

<div id='id-tips'/>

---

<div id='id-cleanr'/>

### CLEAN-R: testing intermodal associations

```R
mod=model.matrix(~covariates)
fit = CleanR(ymat = data1, xmat = data2, mod = mod, distmat = distmat)
```

Please refer to the manual for more information about the optional arguments.

```R
help(CleanR)
```


---

<div id='id-cleanv'/>

### CLEAN-V: testing reliability or heritability (forthcoming)

(Forthcoming)

---

<div id='id-cleanvisualize'/>

### Visualization 

(Forthcoming)


### Frequently asked questions:
<div id='id-q1'/>

**How do I extract surface data from HCP?**
Please refer [ciftiTools](https://github.com/mandymejia/ciftiTools). Once you install [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) in your computer and obtain data files in cifti format and surface information in the `surf.gii` format, then 

```R
library(ciftiTools)
library(rgl)
ciftiTools.setOption("wb_path", "/Applications/workbench")
xii = read_cifti(fname, surfL, surfR, resamp_res = 10242)  #resamp_res: how many vertices to resample
```

Then you can access the cortical data by the following.
```R
xii$data$cortex_left
xii$data$cortex_right
```
The corresponding mesh information can be assessed by the following.
```R
xii$surf$cortex_left
xii$surf$cortex_right
```

<div id='id-q2'/>

**Which surface should we use for registration?**

We use spherical surface as a default and use inflated or midthickness surface for visualization ([reference](https://doi.org/10.1016/j.neuroimage.2016.05.038)). However, please note that there is no definitive answer for this, and there are recent [articles](https://doi.org/10.1016/j.neuroimage.2022.118908) that supported the midthickness surface for registration. Whenever possible, we recommend to conduct an exploratory analysis to make sure the parametric kernel agrees with empirical data. 

<div id='id-q3'/>

**How do I obtain a pairwise distance matrix?**

We recommend using geodesic distance for mesh surfaces. You may use [Python](https://pypi.org/project/pygeodesic/) or [C++](https://code.google.com/archive/p/geodesic/wikis/ExactGeodesic.wiki) to obtain a pairwise geodesic distance matrix.

It requires the extraction of `vertices` and `faces` matrices. The `vertices` matrix contains the 3D coordinate for each vertex, and each row of the `faces` matrix contains the indices of three vertices that construct a triangle in the mesh surface.

If you use [ciftiTools](https://github.com/mandymejia/ciftiTools), it can be accessed directly from the surface GIFTI. If you use [freesurferformats](https://cran.r-project.org/web/packages/freesurferformats/index.html), the surface file (e.g. `lh.pial` or `lh.inflated` from FreeSurfer) can be loaded directly using the `read.fs.surface()` function. Once you load `vertices` and `faces', the following manipulation is necessary due to the difference between R and Python.

```R
surf$faces = surf$faces-1
```

Once you loaded `vertices` and `faces` in Python, the following would provide a pairwise distance matrix.

```python
num_vts = len(vertices)
distmat = np.zeros((num_vts, num_vts))
target_indices = np.array(range(num_vts))
for i in range(num_vts):
    dists, best_source = geoalg.geodesicDistances(np.array([i]), None)
    distmat[i, :] = dists
```

The last step is to subset the distance matrices with your interest, for example, by excluding non-cortex vertices (e.g. a medial wall), which should be straightforward. 

<div id='id-q4'/>

**Is it possible to fit CLEAN/CLEAN-R/CLEAN-V separately for two hemispheres and combine results afterwards?**

Yes, it is **necessary** to set a brain-wise threshold that controls FWER at the nominal level. Please make sure you specify the same seed (`seed`) and the same number of resamples (`nperm`) for the `Clean()` function. Then you may use the `combine()` function to get a new threshold.

```R
Clean.fit.lh = Clean(dataLH, distmatLH, nperm = 5000, seed = 1)
Clean.fit.rh = Clean(dataRH, distmatRH, nperm = 5000, seed = 1)
Clean.fit.combine = combine(Clean.fit.brain = list(Clean.fit.lh, Clean.fit.rh), alpha = 0.05)
```

<div id='id-q5'/>

**What is the recommended value for** `max.radius`**?**

The `max.radius` determines the degree of the spatial domain you're borrowing information from. Higher sensitivity obtained from a large value of `max.radius`, however, comes with the cost of decreased specificity. It should be determined a priori before obtaining any result. We empirically found values between 10 and 20 useful for interpretation. Please report these values in your article.
  
 <div id='id-misc'>

---

### Miscellaneous
Please check out the [SpLoc](https://github.com/junjypark/SpLoc) package, a close family of CLEAN, that conducts clusterwise inference for longitudinal neuroimaging data in comparing two groups's growth/decay. It currently does not support leveraging spatial autocorrelations.

> Park JY, Fiecas M (2021) Permutation-based inference for spatially localized signals in longitudinal MRI data. Neuroimage, 239, 118312. [article link](https://doi.org/10.1016/j.neuroimage.2021.118312)


