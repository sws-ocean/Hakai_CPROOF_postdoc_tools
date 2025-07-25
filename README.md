# Oceanographic Tools: Northeast Pacific Observational Analyses

This repository provides MATLAB and python scripts and resources developed primarily for analyzing observational oceanographic data in the Canadian Pacific, with a focus on C-PROOF glider missions, DFO CTD and mooring datasets, Hakai Institute CTD datasets, and water mass analyses. These tools were developed during research projects conducted across Queen Charlotte Sound, the Central Coast, and the broader northeast Pacific. Some scripts are a little untidy, and several contain legacy code blocks near the end of the file, which are retained for  context but may not be actively used.

## 🔍 Contents

Each subdirectory contains scripts, metadata, and documentation specific to a project:

- [`glider_processing/`](glider_processing/):  
  Tools for working with regional ocean glider data, including preprocessing, oxygen sensor corrections, and structuring for analysis.

- [`dfo_processing/`](dfo_processing/):  
  Scripts for accessing and processing CTD, bottle, and mooring data from Fisheries and Oceans Canada (DFO) databases.

- [`qcs_oxygen/`](qcs_oxygen/):  
  Analysis code for the manuscript "Oxygen variability on the Canadian Pacific shelf: trends, drivers, and projections in the context of emerging hypoxia in Queen Charlotte Sound".

- [`centralcoast/`](centralcoast/):  
  MATLAB implementation of Optimum Multiparameter (OMP) analysis for attributing water mass contributions in Central Coast datasets.

- [`np_circulation/`](np_circulation/):  
  Scripts for analyzing large-scale circulation patterns in the northeast Pacific related to the Canadian Pacific.

## 🛠 Requirements

- MATLAB Toolboxes:
  - Gibbs SeaWater Toolbox
  - m_map toolbox

## 📂 Data Access

This repository does **not** include raw data. Where possible, instructions for accessing the relevant datasets (e.g., DFO, Hakai, C-PROOF) are provided in each subfolder’s README.


## 👤 Author

This repository was developed by Sam Stevens during a postdoctoral fellowship at the Hakai Institute.
