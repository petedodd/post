#!/bin/bash
# script to run all analyses
# see details in the 000-file for expected file structure

# lifetable work
R --slave < 000_lifetables.R

# 00_utilities.R are utilities used by the below scripts

# analysis of treated cases
R --slave < 01_interppop.R
R --slave < 02_HIVinTB.R
R --slave < 03_TBtxOutcomes.R
R --slave < 04_TBnotes.R
R --slave < 05_CDRprep.R
R --slave < 06_agedisag.R
R --slave < 07_notecombine.R
R --slave < 08_treated.R

# analysis of untreated cases
R --slave < 09_incidence.R
R --slave < 10_CFRs.R
R --slave < 11_untreated.R

# make figures and tables
R --slave < Figures.R
R --slave < Tables.R

