netcdf
======

Interpolate dated polygon map series to a regular 3-d grid (x,y,t) stored as a netcdf variable

1. Create NetCDF file with x,y,t dimensions
2. Sample dated polygons at grid coordinates to create new x,y,t variable
3. sample raster at x,y grid coordinates to create new x,y variable

#### Usage

`ncgrid.py [-h] {Create,Sample,Raster} ...`

* **Create** - Create a NetCDF grid from an input raster and date sequence.
    * **ncfile** - output netcdf path
    * **raster** - input raster path
    * **start** - start date YYYY-MM-DD
    * **end** - end date YYYY-MM-DD
    * **count** - count of date coordinates to be spaced equally from start to stop inclusive

* **Sample** - Sample dated polygon map series to netcdf grid using a query. For each date in the time dimension select polygons where *polyDate* is *ineq* *gridDate* extract *targetfield* choosing max *priorityfield* if n>1')
    * **ncfile** - netcdf file path
    * **varname** - created variable name
    * **infeatures** - input polygon shapefile path
    * **target** - extracted variable field name
    * **ineq** - inequality used in query
    * **priority** - priority field if multiple polygons meet query
    * **datefield** - input date field name
    * **dateformat** - input date format string (strptime format)
    * **dtype** - created variable data type
    * **revpri** - (optional) reverse order indicated by priority field (-1 * x)

* **Raster** - Add a single raster variable to existing netcdf file
    * **ncfile** - netcdf file path
    * **varname** - created variable name
    * **raster** - input raster
    * **dtype** - data type
