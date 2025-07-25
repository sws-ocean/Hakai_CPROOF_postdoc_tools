# Central Coast Oxygen Analysis and Seasonal Trends

This folder contains scripts and figures related to oxygen variability and water mass analysis in the Central Coast of British Columbia, with a focus on the deep renewal layer in Fitz Hugh Sound and Burke, Dean, and Rivers Inlets.

The code includes tools for solving the Optimum Multiparameter (OMP) mixing problem, visualizing seasonal oxygen patterns, and comparing deep water renewal characteristics across the Central Coast region.

## Contents

| File                    | Description |
|-------------------------|-------------|
| `OMP_FHZ08.m`           | Performs OMP analysis for station FHZ08 mooring. Solves for water mass fractions and residuals. |
| `fit_harmonics.m`       | Fits annual harmonics to time series data (e.g., oxygen) using least squares. Used in seasonal analysis. |
| `sinfit_Calvert.m`      | Performs sinusoidal fits to seasonal oxygen cycles across multiple stations. Generates seasonal parameter maps. |
| `CCDWR.m`               | Computes Central Coast Deep Water Renewal metrics from temperature, salinity, and oxygen data. |
| `CCDWR.jpg`             | Schematic showing water mass structure during renewal events. |
| `CCDWR_TnO.jpg`         | Plot showing temperature and oxygen anomalies during renewal events. |
| `seasonalityCalvertAnnotate.png` | Annotated map of seasonal oxygen variability (subset). |
| `seasonalityCalvertBnDAnnotate.png` | Annotated map of seasonal oxygen variability across all stations. |

## Description of Analysis

- **OMP Analysis** (`OMP_FHZ08.m`) estimates contributions of source water masses to observations at the FHZ08 station using salinity, temperature, and oxygen as tracers. Results help identify seasonal variability in source fractions.
  
- **Seasonal Analysis** (`sinfit_Calvert.m`) quantifies the amplitude and phase of seasonal oxygen cycles across a network of moorings using sinusoidal curve fitting. Plots illustrate spatial differences in timing and strength of subsurface oxygen renwal.

- **CCDWR Metrics** (`CCDWR.m`) defines the timing and extent of Central Coast deep water renewal events. 

## Dataset
- When using Hakai Institute CTD data, I tend to download .mat files from the [Hakai Institute ERDDAP server](https://catalogue.hakai.org/erddap/index.html) and work from them. DFO CTD data are sourced from the IOS [waterproperties.ca](waterproperties.ca) archive [see this README for DFO data workflow](../dfo_processing/README.md)

## Requirements

- [OMP package by Johannes Karstensen and Matthias Tomczak](https://omp.geomar.de/README.html)

---
