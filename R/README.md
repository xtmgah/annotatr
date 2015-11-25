## Description of `R/`

This folder contains all `R` code/functions that the user will have access to. It can also contain helper R code/functions. It would be most helpful to split the functions among files by what the function does. For example, functions that read user data could go in `read.R`, functions that visualize data could go in `visualize.R`, etc.

Unit tests for the code/functions contained in this folder belong in `tests/testthat`. Tests should likewise be split among files by what function is being tested.
