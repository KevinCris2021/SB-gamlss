# SB-gamlss
Simplex Binomial (SB) family for GAMLSS in R

## Quick start

1) Download this repository (Code -> Download ZIP) or clone it.
2) In R, set your working directory to the repository root folder.
3) Install the package from GitHub and load it in R:

```r
# install.packages("remotes")
remotes::install_github("KevinCris2021/SB-gamlss")

library(SBgamlss)
library(gamlss)

fit <- gamlss(
  y ~ x,
  sigma.fo = ~ 1,
  family = SB(),
  data = df,
  bd = df$bd
)

