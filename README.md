netcdf
======

Interpolate dated polygon map series to a regular 3-d grid (x,y,t) stored as a netcdf variable


## Usage

`ncgrid.py [-h] {Create,Sample,Raster} ...`

Command | Description           |
--------|-----------------------|
Create  | Create NetCDF file    |
Sample  | Sample dated polygons |
Raster  | sample xy raster      |


### **Create**

Create a NetCDF grid from an input raster and date sequence.

* **ncfile** - output netcdf path
* **raster** - input raster path
* **start** - start date YYYY-MM-DD
* **end** - end date YYYY-MM-DD
* **count** - count of date coordinates to be spaced equally from start to stop inclusive

### **Sample**

Sample dated polygon map series to netcdf grid using a query. For each date in the time dimension select polygons where *polyDate* is *ineq* *gridDate* extract *targetfield* choosing max *priorityfield* if n>1')

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

### **Raster**

Add a single raster variable to existing netcdf file

* **ncfile** - netcdf file path
* **varname** - created variable name
* **raster** - input raster
* **dtype** - data type

### Output file

* dimensions x,y,t
* coordinate variables xcoord, ycoord, tcoord, datestr 
* additional attributes
    * spatial_reference
    * linear_unit
    * cellsize

