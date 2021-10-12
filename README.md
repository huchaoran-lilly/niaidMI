# README #

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/niaidMI)](https://cran.r-project.org/package=niaidMI)

This README would normally document whatever steps are necessary to get this application up and running.

### What is this R package for? ###

* This package provides the implementation of **Multiple Imputation with Markov Chain** in application to a popular COVID-19 scale, **NIAID-OS**.

* Version: 1.0.0

### How do I get set up? ###

* Run the following R code.

```
if (!require("devtools")) install.packages("devtools")
devtools::install_github("huchaoran-lilly/niaidMI")
```

Note: This package includes Cpp codes, so you will need a CPP compiler to install.
You can install necessary tools by using Rtools provided by [CRAN](https://cran.r-project.org/).

* Run `library(niaidMI)` to load in R.

### How do I use? ###

Please check the documentation and help file within this package.

### Who do I talk to? ###

* Nathan Morris, <morris_nathan@lilly.com>
* Chaoran Hu, <hu_chaoran@lilly.com>

