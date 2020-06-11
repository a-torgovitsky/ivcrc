# ivcrc: An Instrumental Variables Estimator for the Correlated Random Coefficients Model

**Description**: A Stata module for implementing the instrumental variables correlated random coefficients estimator proposed in Masten and Torgovitsky (2016) and analyzed in more detail in Masten and Torgovitsky (2014).  The module also allows users to estimate varying coefficient models.

**Authors**: This module was written by [David Benson](https://www.federalreserve.gov/econres/david-a-benson.htm), in collaboration with [Matt Masten](http://www.mattmasten.com) and [Alexander Torgovitsky](https://sites.google.com/site/atorgovitsky/).

## Requirements

* Stata version 12

## Installation

In Stata, type "ssc install ivcrc".  Or copy `ivcrc.ado`, `_ivcrc_estimator.ado`, and `ivcrc.sthlp` files to your Stata/ado/personal directory. See the Stata help file for explanation of the syntax. 

## Subdirectories

* Application - Files for our empirical illustration
* Test - Files for testing correctness of the module's output

## Troubleshooting

Please post problems or suggestions to the issue queue.

## References

Masten and Torgovitsky (2016) [Identification of Instrumental Variable Correlated Random Coefficients Models](https://doi.org/10.1162/REST_a_00603), _The Review of Economics and Statistics_ 98 (5), pp. 1001-1005

Masten and Torgovitsky (2014) [Instrumental Variables Estimation of a Generalized Correlated Random Coefficients Model](http://www.cemmap.ac.uk/wps/cwp021414.pdf), Cemmap working paper CWP02/14

## License

&copy; 2018 David Benson, Matt Masten, Alexander Torgovitsky

The contents of this repository are distributed under the MIT license. See file
`LICENSE` for details.

