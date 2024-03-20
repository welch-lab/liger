# rliger 2.0.0

## R CMD check results

There were no ERRORs or WARNINGs.

There were two NOTEs from winbuilder:

```
* checking CRAN incoming feasibility ... [86s] NOTE
Maintainer: 'Yichen Wang <wayichen@umich.edu>'

Suggests or Enhances not in mainstream repositories:
  DoubletFinder, RcppPlanc
Availability using Additional_repositories specification:
  DoubletFinder   yes   https://blaserlab.r-universe.dev
  RcppPlanc       yes   https:/welch-lab.r-universe.dev 
```

```
* checking package dependencies ... NOTE
Packages suggested but not available for checking:
  'DoubletFinder', 'RcppPlanc'
```

*DoubletFinder* and *RcppPlanc* are hosted on R-universe and used conditionally in rliger. 

Besides, we kindly request that CRAN team can approve our submission for *RcppPlanc* so we can move it to Imports.
