# iHPlots

Created by Ruth Gómez Graciani

## About

Give project objective, description and scope

## Structure

* `analysis`: a directory for each new run should be created, with the parameter files for that run if necessary.
* `code`: `software/` is for external code and the others are original code.
* `data`: `raw/` for data as-is and `use/` for preprocessed data that is the same regardless of analysis parameters. Please explain origin, content and format of each file in `raw/` to make it reproducible. 
* `project`: papers and other useful bibliography, project history.
* `report`: reports, power points, etc. each in its directory. Final, sent files should be in deliver/ subdirectories.
* `tmp`: temporary files, organized as in `analysis`. 

## Prerequisites

This project should be associated with a conda environment, which can be imported from the .yml file in the main directory.

## Running the tests


The analysis is divided in different parts, that are options in `runall.sh`. Each of them could have its own requirements.

## A note on code comment structure

The header: 
```
#!/bin/bash
# Ruth Gómez Graciani
# 16 11 2020

###############################################################################
# Description:                                                                
#                                         
###############################################################################
```

The levels of code:

```
# MAIN BLOCK
# Description
# =========================================================================== #

# MAIN BLOCK - SUB BLOCK
# --------------------------------------------------------------------------- #

# Normal comment

```

Note the 4-3-2-1 lines in the different comments.

