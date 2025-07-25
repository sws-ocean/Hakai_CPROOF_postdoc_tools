# DFO Data Access and Processing Tools

This folder contains MATLAB scripts for downloading, loading, and structuring CTD and mooring data from Fisheries and Oceans Canada (DFO), particularly from the [IOS Water Properties archive](https://www.waterproperties.ca/). These tools support analysis of hydrographic stations and time series from moored instruments (e.g., ADCPs and CTDs).

Some of the scripts here are adapted from scripts written of Rich Pawlowicz and Ben O’Connor-thanks to them both. 

## Contents

| File                 | Description |
|----------------------|-------------|
| `ios_rd.m`           | Reads IOS format files into MATLAB arrays. This routine was written by Rich Pawlowicz. |
| `ios_ctd.m`          | An example of how to process groups of IoS format files.
| `loadMooringADCP.m`  | Loads and reformats ADCP time series data from DFO moorings. Converts to time-depth grids with QC steps. |
| `loadMooringCTD.m`   | Loads individual Seabird or RBR moored CTD files from a directory, applies standard QC and depth referencing. |
| `readMooringCTDData.m` | Reads metadata and time series from DFO moored CTD deployments. Includes time parsing and pressure referencing logic. |

---
