# spatTempBFA
This repository implements novel practicable spatiotemporal Bayesian Gaussian factor models detailed in the manuscript "Enhancing Scalability in Bayesian Nonparametric Factor Analysis of Spatiotemporal Data" by Cheng and Li (2023) `<arXiv:2312.05802>`. The baseline framework is a separable factor model with temporally dependent factors and spatially dependent factor loadings from a spatial mixture model under Probit Stick Breaking Processes (PSBP), which enables spatial clustering of temporal trends. The covariance structure on the latent temporal factors can be specified as the AR(1) process `ar1`, the exponential process `exponential`, their seasonal counterparts `sar1`, `sexponential`, or the VAR(1) process (not in the aforementioned manuscript but will be covered in future research; not completely implemented yet in this current preliminary version of the package). We incorporate a new slice sampling algorithm that allows unknown varying numbers of spatial mixture components across all factors and leads to significantly enhanced computational and storage efficiency and scalability. We further integrate into our spatiotemporal models a nearest-neighbor Gaussian process (NNGP) coupled with novel sequential updating schemes for the spatially varying latent variables in the PSBP prior, thus achieving high spatial scalability. Crucial inferential procedures including out-of-sample predictions at future time points or new spatial locations and clustering spatial locations into regions with similar temporal trajectories are supported by our package. A more comprehensive version of the package, which will enable non-homogeneous temporal structure generalizations for the latent factors, is forthcoming.









