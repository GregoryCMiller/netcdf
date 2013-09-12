netcdf
======


## Usage

`ncgrid.py {Create,Sample,Raster} ...`


Subcommand | Description                                                           |
---------- | --------------------------------------------------------------------- |
Create     | Create regular x,y,t netcdf file from an input raster and date series |
Sample     | Sample dated polygon map series at x,y,t grid coordinates.            |
Raster     | Sample x,y raster at netcdf x,y grid coordinates                      |


## **Create**

Create regular x,y,t netcdf file from an input raster and date series.

Parameter | Description                                      |
--------- | ------------------------------------------------ |
ncfile    | output netcdf path                               |
raster    | input raster path                                |
start     | start date YYYY-MM-DD                            |
end       | end date YYYY-MM-DD                              |
freq      | sampling frequency [DAILY,WEEKLY,MONTHLY,YEARLY] |

## **Sample**

Sample dated polygon map series to netcdf grid using a query. 

Parameter  | Description                                     |
---------- | ----------------------------------------------- |
ncfile     | netcdf file path                                |
varname    | created variable name                           |
infeatures | input polygon shapefile path                    |
target     | extracted variable field name                   |
ineq       | inequality used in query                        |
priority   | priority field if multiple polygons match query |
datefield  | input date field name                           |
dateformat | input date format string (strptime format)      |
dtype      | created variable data type                      |
--revpri   | reverse priority (-1 * priority)                |

Query 

    For each date in the time dimension
    select polygons where *polyDate* is *ineq* *gridDate*
    extract *targetfield* (choosing max *priorityfield* if n>1')
    
## **Raster**

Add a single raster variable to existing netcdf file

Parameter| Description           |
-------- | --------------------- |
ncfile   | netcdf file path      | 
varname  | created variable name |
raster   | input raster          |
dtype    | data type             |
