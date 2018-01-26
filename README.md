# networkR

Development version of the R package **networkR**. This package
contains a collection of various semi-useful functions for network af family analysis.

To install the development version of **networkR** run the following
command from within R (this requires that the devtools package is
already installed on the system.)

```r
devtools::install_github('ekstroem/networkR')
```

(CRAN checks : [![Travis-CI Build Status](https://travis-ci.org/ekstroem/networkR.svg?branch=master)](https://travis-ci.org/ekstroem/networkR))


## Package overview


### Families

*  `cluster_families` will identify clusters (families) based on information of ids, father ids and mother ids.


### Network analysis


*  `adjacency` 
*  `hits`

### Formats

The `linkage` format for family data ...



The `trio` format that we use in this package is a reduced/restricted version
of the linkage format consisting of three *integer* (numeric) values for each individual: 

*   `id` - the id of the individual,
*   `fid` - the id of the father,
*   `mid` - the id of the mother.

For the `fid` and `mid` variables a 0 or `NA` represents no
information, and a person with nor information on father and mother id
is considered a founder. For consistency, an individual should not
have just a single parent missing.

Note that since we do not have a family indicator all id's should be uniques across all families.

