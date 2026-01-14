# SB-gamlss
Simplex Binomial (SB) family for GAMLSS in R

## Quick start

1) Download this repository (Code -> Download ZIP) or clone it.
2) In R, set your working directory to the repository root folder.
3) Run:

```r
source("R/load_all.R")

library(gamlss)
fit <- gamlss(y ~ x, sigma.fo=~1, family = SB(), data = df, bd = df$bd)
