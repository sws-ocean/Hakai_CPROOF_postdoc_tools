# Oceanographic Tools: Northeast Pacific Observational Analyses

This repository provides MATLAB scripts and resources developed primarily for analyzing observational oceanographic data in the Canadian Pacific, with a focus on C-PROOF glider missions, DFO CTD and mooring datasets, Hakai Institute CTD datasets, and water mass analyses. These tools were developed during research projects conducted across Queen Charlotte Sound, the Central Coast, and the broader northeast Pacific. Some scripts are a little untidy, and several contain legacy code blocks near the end of the file, which are retained for  context but may not be actively used.

## 🔍 Contents

Each subdirectory contains scripts, metadata, and documentation specific to a project:

- [`glider_processing/`](glider_processing/):  
  Tools for preprocessing regional ocean glider data, including oxygen sensor corrections and data structuring.

- [`dfo_processing/`](dfo_processing/):  
  Scripts for processing CTD, bottle, and mooring data from Fisheries and Oceans Canada (DFO) databases.

- [`qcs_oxygen/`](qcs_oxygen/):  
  Analysis code for the manuscript "Oxygen variability on the Canadian Pacific shelf: trends, drivers, and projections in the context of emerging hypoxia in Queen Charlotte Sound".

- [`centralcoast/`](centralcoast/):  
  Scripts and figures related to oxygen variability and water mass analysis in the Central Coast of British Columbia, with a focus on the deep renewal layer in Fitz Hugh Sound and Burke, Dean, and Rivers Inlets.

- [`np_circulation/`](np_circulation/):  
  Scripts and PDF reports analyzing large-scale circulation in the northeast Pacific, with emphasis on water mass pathways influencing the Canadian Pacific.

## 🛠 Requirements

- I frequently use these MATLAB toolboxes:
  - Gibbs SeaWater Toolbox
  - m_map toolbox
  - cmocean toolbox

## 📂 Data Access

This repository does **not** include raw data. Where possible, instructions for accessing the relevant datasets (e.g., DFO, Hakai, C-PROOF) are provided in each subfolder’s README.


## 👤 Author

This repository was developed by Sam Stevens during a postdoctoral fellowship at the Hakai Institute.
