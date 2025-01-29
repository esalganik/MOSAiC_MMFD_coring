A forcing dataset for single column model based on MOSAiC first-year ice (FYI) and second-year ice (SYI) coring data

The dataset covers observations from 29 Oct 2019 to 29 July 2020.
Vertical spatial resolution is 0.05 m, and temporal resolution is 6 h.
The dataset provides sea-ice salinity, temperature, brine volume, density, and snow thickness data.

FYI temperature was measured by ice mass balance buoy 2019T66: https://doi.org/10.1594/PANGAEA.938134
FYI salinity and density were measured by ice coring: https://doi.org/10.1594/PANGAEA.971385
SYI temperature was measured by ice mass balance buoy 2019T62: https://doi.org/10.1594/PANGAEA.938134
SYI salinity and density were measured by ice coring: https://doi.org/10.1594/PANGAEA.959830
Pre-melt snow thickness was measured by various ice mass balance buoys: https://doi.pangaea.de/10.1594/PANGAEA.973193

Files:
T66.mat, T62.mat - ice mass balance buoys raw data
Coring_old_fb.mat - ice coring raw salinity data
SIMBA/2019T**_icethick.tab - ice mass balance buoys processed from Preußer et al., 2025
single_column.m - raw data processing and export to NetCDF
MOSAiC_FYI_coring.nc, MOSAiC_SYI_coring.nc - processed NetCDF files
