[![gammapy](https://img.shields.io/badge/powered%20by-gammapy-orange.svg?style=flat)](https://gammapy.org/)

# roi_analysis_gammapy_catalogs
A Python code to search for possible γ-ray counterparts to the target source and to perform the spectral model fitting. This code selects the sources (in the Gammapy source catalogs) within the region of interest (centered in the position of the target source) and finds the best fit for the given spectrum model.

#### The set_analysis.dat file accepts the following options:

source_name: Source name based on J2000 coordinates  
pos_ra: Right ascension (in degrees) 
pos_dec: Declination (in degrees)
radius_roi: The maximum angle (in degrees) of separation between the target source and its counterparts 
e_ref_min: The minimum reference energy (in TeV) of the flux points table (default is None)
e_ref_max: The maximum reference energy (in TeV) of the flux points table (default is None)

