# ivcrc
A Stata module for implementing the instrumental variables correlated random coefficients estimator proposed in Masten and Torgovitsky (2016, published in *The Review of Economics and Statistics* 98 (5), pp. 1001--1005) and analyzed in more detail in Masten and Torgovitsky (2014, Cemmap working paper CWP02/14).

The module was written by David Benson (Northwestern University).

A help file has been added to assist with the syntax and estimation options. 

The module has two ado files which work together: IVCRC.ado is a wrapper, feeding commands to IVCRC_ESTIMATE.ado which performs the computations. Users must load both for the module to work, either by running them manually or by placing the files in the local ado directory. There are two available versions of the estimation routine IVCRC_ESTIMATE, users choosing to load the files into a local ado directory must rename their preferred routine to IVCRC_ESTIMATE.ado. IVCRC_ESTIMATE_ONEVAR.ado is for the case of a single basic endogenous variable and is faster. IVCRC_ESTIMATE_MULTIVAR.ado is for the case of multiple basic endogenous variables. Both allow arbitrary "derived endogenous" variables.  


