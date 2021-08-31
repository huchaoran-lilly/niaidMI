# OSS Project Description Template

## OSS Project Title

R package `niaidMI`

## Contributors

* Nathan Morris, morris_nathan@lilly.com
* Chaoran Hu, hu_chaoran@lilly.com

## Local Business Area

* Jingyi Liu, Sr Director-Biostatistics, Statistics-Immunology, liu_jingyi@lilly.com

* Brenda J Crowe, Sr Research Advisor, Statistics-Immunology, crowe_brenda_j@lilly.com

## License

The MIT License

```
Copyright 2021 Eli Lilly and Company

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

## Repository

* Development repository: Chaoran Hu's Lilly GitHub: https://github.com/huchaoran-lilly/niaidMI

* Target publish place: CRAN


## Purpose

This package is the implementation of Markov Model Multiple Imputation (MMMI)
with the application to NIAID-OS, a common COVID-19 scale. The daily scores
of NIAID-OS are modeled by a Markov chain and the imputation is processed
based on the fitted model. MMMI incorporates Markov chain to multiple imputation
framework. The application of MMMI can be extended to other discrete scales.

## Implementation

`niaidMI` is implemented as a typical R package, with C++ code to improve
its performance. It depends on two external R packages, `Rcpp` and `testthat`.

## Development

* Development and maintenance of the code will continue for the foreseeable future.

* Code development and maintenance will be the responsibility of Nathan Morris and Chaoran Hu.

* Automated checks are and will continue to be used as the project is updated.
  Independent validation has been performed to assure the high quality of code.

* External contributions (e.g. pull requests) will be evaluated before merging them into the project.

