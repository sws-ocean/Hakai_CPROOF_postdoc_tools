# Glider Data Processing and Oxygen Correction

This folder contains MATLAB code and calibration files used to process delayed-mode glider observations from the Northeast Pacific. It includes tools for loading, organizing, and correcting glider data—particularly dissolved oxygen—from multiple sensor types (e.g., AADI, AROD-FT, JFE).

## Contents

| File/Folders      | Description |
|-------------------|-------------|
| `Gdataset.m`      | Loads and standardizes gridded delayed-mode glider data into a structured dataset (`allProfs`), including georeferencing onto a canyon axis. |
| `GoxyCore.m`      | Main script for performing multi-stage oxygen correction (P, T, S, and response-time correction) for AADI, AROD, and JFE sensors. Also compares corrected glider oxygen with nearby CTD profiles and applies drift correction. |
| `Gstats.m`        | Quick diagnostic summary of the glider dataset (e.g., number of missions, mean profile spacing/frequency, histograms). |
| `rdAROD.m`        | Reads oxygen sensor calibration data from `AROD_cal.csv`. |
| `rdAADI.m`        | Reads oxygen sensor calibration data from `AADI_cal.csv`. |
| `AROD_cal.csv`    | Calibration coefficients for AROD-FT oxygen sensors. |
| `AADI_cal.csv`    | Calibration coefficients for AADI oxygen sensors. |

## Workflow

1. **Build Dataset**  
   Use `Gdataset.m` to aggregate and standardize glider `.nc` files into the `allProfs` structure. This step also:
   - Converts time formats
   - Maps data onto a defined canyon transect
   - Computes TEOS-10 variables (SA, CT, N²)

2. **Apply Oxygen Corrections**  
   Use `GoxyCore.m` to apply:
   - Pressure and salinity corrections
   - Response-time (tau) corrections using phase lag analysis
   - Optional drift correction via regression with CTD profiles
   NOTE: this is provided as a reference, not as a definitive guide on oxygen sensor corrections. I referenced Hayley Dosser's [processing reports](https://cproof.uvic.ca/gliderdata/deployments/reports/) to write this code.

3. **Quality Control and Visualization**  
   - `Gstats.m` provides basic statistics and spacing/frequency plots
   - Output includes diagnostic figures for each mission
   - Corrected oxygen is stored in `allProfs.oxygen_corrected`

## Input Requirements

- Gridded delayed-mode `.nc` C-PROOF glider data available [on the C-PROOF website](https://cproof.uvic.ca/deployments/index.html)
- Calibration CSVs for each oxygen sensor type
- TEOS-10 and plotting packages:
  - `gsw` (Gibbs Seawater Toolbox)
  - `m_map`
  - `export_fig` (optional)

## Outputs

- `glider_allProfs.mat`: main structured dataset used in further analysis
- Corrected oxygen field: `allProfs.oxygen_corrected`
- Figures for QA/QC of each mission in `quickCorrFigs/` (user must create this directory if saving locally)
- Optionally: `corrected_allProfs.mat` and `Dcorrected_allProfs.mat` for intermediate stages

## Notes

- Sensor-specific correction routines depend on matching serial numbers to the appropriate calibration sheet.
- Code is configured for C-PROOF Calvert line deployments; paths may need to be updated for different systems.
- Time-domain oxygen sensor corrections are based on the tau-corrected phase method adapted from Uchida et al. (2010). See Hayley Dosser's [processing report](https://cproof.uvic.ca/gliderdata/deployments/reports/) for reference.


---

