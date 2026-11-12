# ------------------------------------------------------------------------------------------
# Determining the Best DNA Reference Database for Zooplankton in the San Francisco Estuary
# 
# Kristina Nguyen
# Nov 2026
# ------------------------------------------------------------------------------------------

# load libraries and functions
library(tidyverse)
library(data.table)

source("compare-reference-databases/compare-two-databases.R")
source("compare-reference-databases/compare-three-databases.R")

# import data (BLAST results, contain all columns)
midori <- fread("", fill = T, col.names = c("query_seqid", "result_seqid", "evalue", "bitscore", "length",
                                            "percent_coverage", "nident", "percent_identity"))

# standardize taxonomy

# We also need to add higher taxonomy for the resulting matches. Our functions require the columns 
# kingdom, phylum, class, order, family, genus, and species. For databases, like PR2, which use 
# non-traditional taxonomic levels, we standardize the taxnomic names using a more traditional database, 
# like GBIF, before running these functions. The optional column, sciname, is the taxonomic name of the 
# lowest resolution. For example, if a match was resolved only to the family level, sciname would have 
# the family name, while genus and species would have NA values. If sciname is not included in the inital 
# dataframe, it will be added through addScinameLevel() within the compare functions (compareBOLDandMIDORI, 
#                                                                                     compareMLMLandLavrador, etc.).

# keep best hit per database using best shared method (Pappalardo 2025)

# best-hit: we kept the top hit from the BLAST results, regardless of whether the second hit had identical 
# similarity measures.
# best-shared: when two or more rows had identical similarity measures (evalue, bitscore, length, qcovs, 
#                                                                       nident and pident), we used a customized 
# function to only keep the taxonomic assignment (or taxonomic level) for which the matches agreed. For example, 
# identical matches to Balanus glandula and Balanus amphitrite will only end up with classification to Balanus; 
# but identical matches to Balanus glandula will retain the species Balanus glandula.

# comparisons

bold_vs_midori <- compareBOLDandMIDORI(your_midori_df, your_bold_df, identiy_th) %>%
  pickFinalTax_BoldVsMidori()

pr2_vs_silva <- compareBOLDandMIDORI(your_pr2_df, your_silva_df, identiy_th) %>%
  pickFinalTax_BoldVsMidori()

mlml_vs_lavrador <- compareMLMLandLavrador(your_mlml_df, nagarajan, identiy_th) %>%
  pickFinalTax_MLMLvsLavrador()

your_final_df <- compare3Databases(_, _, _, idne) %>%
  pick3FinalTax()

# standardize scientific names to single taxonomic framework (GBIF or WoRMS)
