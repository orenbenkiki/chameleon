Chameleon
=========

Overview
--------

`chameleon` contains a set of functions for automatically assigning colors to multi-dimensional
data. The colors are distinct (as much as possible) and also reflect the internal structure of the
data (as much as possible).

See the [published R package](https://CRAN.R-project.org/package=chameleon), the
[rdocumentation](https://www.rdocumentation.org/packages/chameleon), or the [latest github version
documentation](https://tanaylab.github.io/chameleon/index.html) for details. Specifically, the
[vignette](https://tanaylab.github.io/chameleon/articles/examples.html) explains why and
how to use this package.

Installation
------------

To install it, use:

``` r
install.packages('chameleon')
```

Usage
-----

In general, if you have a matrix whose rows represent some elements/observations and whose columns
represent some variables/measurements, then use `chameleon::data_colors` to obtain a color for each
row. You can pass this function the optional `groups` parameter which assigns a group to each row,
to obtain a color for each group instead of a color for each row.

See the [published R package](https://CRAN.R-project.org/package=chameleon), the
[rdocumentation](https://www.rdocumentation.org/packages/chameleon), or the latest github version
documentation's [reference section](https://tanaylab.github.io/chameleon/reference/index.html) for
the list and description of the available functions.
