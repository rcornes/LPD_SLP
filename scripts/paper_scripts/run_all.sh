#!/bin/bash

# --------------------------------------------------------------------------------
# Script to produce all of the plots/tables
# Written by R.Cornes on 31-12-2023
# --------------------------------------------------------------------------------

cd ~/projects/geowind2/scripts

# Get 20CR data
# ./download_20CR_jet_ensemble.sh
# ./download_20CR_jet.sh

# Generate Geowind Percentiles
Rscript make_geowind.R
Rscript make_percentiles.R
Rscript make_boot.R noon

# Main Paper plots
Rscript fig_1.R
Rscript fig_2.R
Rscript fig_3.R
Rscript fig_4.R
Rscript fig_5.R

# Supplementary plots/tables
Rscript make_geowind_nonoon.R
Rscript make_percentiles_3hr.R
Rscript make_boot.R 3hr
Rscript make_percentiles_nonoon.R
Rscript geowind_trend_table.R
Rscript fig_2_nonoon.R
Rscript deltap_plots.R
Rscript nao_correlations.R
Rscript noon_vs_3hr_stats.R
Rscript plot_percentiles_annual.R
Rscript plot_annual_frequency.R
Rscript plot_UV.R
Rscript ks_test.R
#Rscript geowind_20CR.R

