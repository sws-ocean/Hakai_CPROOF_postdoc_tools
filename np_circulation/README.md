# North Pacific Circulation

This folder contains scripts and figures related to the influence of large-scale North Pacific circulation on the composition of water masses reaching the Canadian Pacific margin. The focus is on identifying important circulation modes and their relationship to changes in oxygen and water mass structure along the shelf. This work is (yet) unpublished.

## Contents

### Reports

- **`breathing_bifurcation_PCAreport.pdf`**  
  A standalone PDF report documenting an EOF (Empirical Orthogonal Function) analysis of dynamic height anomalies in the North Pacific. It revisits and reproduces Freeland's (2006) breathing and bifurcation modes using dynamic height differences in the N. Pacific.

- **`Can_Pac_water_masses.pdf`**  
  Report describing a mixing analysis on the 26.5 σ₀ isopycnal using temperature–salinity relationships to estimate relative contributions of Pacific Equatorial Water (PEW) and Pacific Subarctic Upper Water (PSUW) at Canadian margin stations. It connects shifts in water mass composition to oxygen variability.
  
### Trajectories
- **`allTraj_26_5.jpg`**
  Backward particle tracking experiments from August 2023 to August 2014 using 3D monthly velocity fields from the ORAS5 ocean reanalysis (1/4° horizontal resolution, 75 vertical levels). Particles are initialized along cross-shelf transects and tracked backward in time using the Parcels Lagrangian framework. Parcels is used to advect particles on the σ₀ = 26.5 isopycnal surface. A stochastic diffusion scheme with horizontal diffusivity set to 100 m²/s is applied to account for unresolved mesoscale motions.

### Scripts

- **`breathing_bifurcation.m`**  
  MATLAB code that reconstructs the NPC (North Pacific Current) breathing and bifurcation modes using dynamic height differences across oceanographic stations.

- **`NEPGindex_ISO.m`**  
   Interpolates reanalysis temperature to the σ₀ = 26.6 isopycnal, performs EOF analysis on resulting anomalies, and compares principal components to regional indices (e.g., spice, bifurcation mode).

- **`subarcticVSequatorial.m`**  
  Projects observed T–S points onto the PEW–PSUW mixing line to calculate fractional contributions from each water mass.

## References

- Freeland, H. J. (2006), *Atmosphere-Ocean*, [https://doi.org/10.3137/ao.440401](https://doi.org/10.3137/ao.440401) 
- Thomson, R. E., & Krassovski, M. V. (2010), *JGR Oceans*, [https://doi.org/10.1029/2010JC006280](https://doi.org/10.1029/2010JC006280)

