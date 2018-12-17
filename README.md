*Product created by Olivier Beauchard, DIVA interpolations by Charles Troupin*

For all the products, the same approach is adopted and described in a jupyter-notebook with the corresponding Julia code:

## Data reading

All data is available on https://dox.ulg.ac.be/index.php/s/EvIwSvTwhtJ6Tmu /benthos_fish_data/.

The text files are read using the `readdlm` function (reading of files where columns are separated by a delimiter, a comma in this case).
In most of the notebooks, the reading functions (as well as the plotting and writing functions) are defined in the corresponding `src` file.

For example, in the notebook `make_fish_product.ipnb`, the functions are defined in `src/make_fishs_products.jl`.

## Interpolation

### Creation of the land-sea mask

The mask is based on the GEBCO bathymetry with a lower resolution than the original, 30-minute resolution product, as such a fine topography is not necessary for the present applications. 

### Parameter choice

Some tests and sensitivity analysis have been performed to find the optimal values of the two main analysis parameters:
1. correlation length and 
2. the signal-to-noise ratio.

As the operations to estimate these parameters are rather costly, we decided to fix these values.

### Data transformation

When the dataset consists of counts of different species, a transformation f(x) = log_10(x + 1.) was applied to the observations. The inverse transformation is then applied to the resulting gridded field obtained in the next step.

### Spatial interpolation

DIVAnd is applied to the transformed data. When necessary, the relative abundance of different species is computed. 

### Error field estimation

The goal of this step is to obtain a field that represents the confidence one can have in the results. In general this field depends on the data coverage, and it can be used to mask the interpolated field when the error is above a given threshold.

## Plotting

Some basic functions were designed for each application, showing the data positions and values (scatter plot) or the interpolated fields, along with the land-sea mask constructed in a previous stage of the analysis.

## NetCDF file writing

The interpolated and error fields are stored in netCDF files.
The name of the file is the same as the data file employed as an input for the interpolation, except that the file extension is changed to ".nc".

Whenever possible, the chosen approach was to write all the information in a single file: for instance, for the fish temporal product (3-year periods for the data), a time dimension was added to the file to that all the periods are available in the same netCDF. For the traits, a dimension `names` and a variable `names` storing the trait names are added.
