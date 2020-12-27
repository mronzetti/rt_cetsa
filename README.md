# RT-CETSA
# Authored by Michael Ronzetti, NIH/NCATS 2020

888888ba  d888888P           a88888b.  88888888b d888888P .d88888b   .d888888
88    `8b    88             d8'   `88  88           88    88.    "' d8'    88
a88aaaa8P'   88             88        a88aaaa       88    `Y88888b. 88aaaaa88a
88   `8b.    88    88888888 88         88           88          `8b 88     88
88     88    88             Y8.   .88  88           88    d8'   .8P 88     88
dP     dP    dP              Y88888P'  88888888P    dP     Y88888P  88     88
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

These r files serve as the beginning of an integrated data analysis package for 
RT-CETSA experimental files. Files are converted first to a format that moltenprot 
(1), an analysis program used forcurve fitting and thermodynamic modeling, accepts. 
After, these returned files are further processed to assignidentifiers and determine
assay statistics.
________________________________________________________________________________________
R Packages
  c('tidyverse',
  'readxl',
  'janitor',
  'Biocmanager',
  'PharmacoGx',
  'ggbeeswarm',
  'cowplot',
  'ggrepel',
  'rlang')
________________________________________________________________________________________
Description of files

1_  label_plate
  |_  Takes in raw matlab files from RT_CETSA experiments and cleans/arranges them for analysis
      in moltenprot software.
2_  process_moltenprot
  |_  Retrieves data from moltenprot analysis folders/files and processes/tidies them for additional analysis.
  
________________________________________________________________________________________
References

(1) Kotov, V, Mlynek, G, Vesper, O, et al. In‐depth interrogation of protein thermal unfolding data with MoltenProt
.Protein Science. 2021; 30: 201– 217. https://doi.org/10.1002/pro.3986