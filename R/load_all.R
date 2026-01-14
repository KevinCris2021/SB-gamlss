# Load and compile SB-gamlss implementation (R + C++)
# Usage:
#   setwd("SB-gamlss")  # repo folder
#   source("R/load_all.R")

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Paquete 'Rcpp' no instalado. Instala con install.packages('Rcpp').")
}

# 1) Compile C++ (creates cpp_joint_sb_ in the R session)
cpp_path <- file.path("src", "SB.cpp")
if (!file.exists(cpp_path)) {
  stop("No se encontró 'src/SB.cpp'. Asegúrate de correr esto desde la raíz del repo.")
}
Rcpp::sourceCpp(cpp_path)

# 2) Source R wrappers and family
source(file.path("R", "dsb_marginal.R"))
source(file.path("R", "dSB.R"))
source(file.path("R", "pSB.R"))

# Optional (if present)
q_path <- file.path("R", "qSB.R")
r_path <- file.path("R", "rSB.R")
if (file.exists(q_path)) source(q_path)
if (file.exists(r_path)) source(r_path)

source(file.path("R", "SB_family.R"))
