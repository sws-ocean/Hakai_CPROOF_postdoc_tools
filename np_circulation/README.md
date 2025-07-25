# North Pacific Circulation

This folder contains scripts and figures related to the influence of large-scale North Pacific circulation on the composition of water masses reaching the Canadian Pacific margin. The focus is on identifying key circulation modes and their relationship to changes in oxygen and water mass structure along the shelf.

## Contents

### 📄 Reports

- **`breathing_bifurcation_PCAreport.pdf`**  
  A standalone PDF report documenting an EOF (Empirical Orthogonal Function) analysis of dynamic height anomalies in the North Pacific. It revisits and reproduces Freeland's (2006) breathing and bifurcation modes using dynamic height differences in the N. Pacific.

- **`Can_Pac_water_masses.pdf`**  
  Report describing a mixing analysis on the 26.5 σ₀ isopycnal using temperature–salinity relationships to estimate relative contributions of Pacific Equatorial Water (PEW) and Pacific Subarctic Upper Water (PSUW) at Canadian margin stations. It connects shifts in water mass composition to oxygen variability.

### 📜 Scripts

- **`breathing_bifurcation.m`**  
  MATLAB code that reconstructs the NPC (North Pacific Current) breathing and bifurcation modes using dynamic height differences across oceanographic stations.

- **`NEPGindex_ISO.m`**  
   Interpolates reanalysis temperature to the σ₀ = 26.6 isopycnal, performs EOF analysis on resulting anomalies, and compares principal components to regional indices (e.g., spice, bifurcation mode).

- **`subarcticVSequatorial.m`**  
  Projects observed T–S points onto the PEW–PSUW mixing line to calculate fractional contributions from each water mass.

### 📘 README

- **`README.md`**  
  This file. Provides documentation for the folder contents and methodology.

## Notes

- All scripts are written in MATLAB.
- GLORYS reanalysis data and observational profiles on isopycnals (e.g., 26.5 σ₀) are assumed as input.
- This folder supports (yet) unpublished work examining links between shelf deoxygenation and offshore circulation variability.

## References

- Freeland, H. J. (2006), *Atmosphere-Ocean*, [https://doi.org/10.3137/ao.440401](https://doi.org/10.3137/ao.440401) 
- Thomson, R. E., & Krassovski, M. V. (2010), *JGR Oceans*, [https://doi.org/10.1029/2010JC006280](https://doi.org/10.1029/2010JC006280)

