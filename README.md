# input4mips-prep
Collection of scripts used to help preparing input4MIPs datasets. These scripts try to manipulate received files as little as possible, hence there is a large reliance on the low-level netCDF4 library and general avoidance of higher-level libraries like xarray and iris which often make other 'helpful' fixes without asking.
