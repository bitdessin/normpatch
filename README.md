[![build](https://github.com/bitdessin/normpatch/actions/workflows/build_pkg.yaml/badge.svg?style=flat-square)](https://github.com/bitdessin/normpatch/actions/workflows/build_pkg.yaml)

# normpatch

*Patch the bias. Keep the biology.*

normpatch is an R package for estimating normalization factors from bulk RNA-seq count data.
It implements methods based on the DEGES framework,
the DEG elimination strategy,
which iteratively reduces the influence of DEGs
when estimating sample-specific scaling factors.
The framework is intended for datasets affected by compositional imbalance
or strongly asymmetric DEG patterns,
where standard normalization assumptions may not hold.



## Documentation

- https://bitdessin.github.io/normpatch/


## Citation

Sun J, Nishiyama T, Shimizu K, Kadota K.
TCC: an R package for comparing tag count data with robust normalization strategies.
BMC Bioinformatics. 14:219, 2013.
[10.1186/1471-2105-14-219](https://doi.org/10.1186/1471-2105-14-219)
