ncgrid
======


## Usage

`ncgrid.py {Create,Sample,Raster} ...`
                                                                                
Subcommand | Description                                                   |
---------- | ------------------------------------------------------------- |
Create     | Create x,y,t netcdf file from an input raster and date series |
Sample     | Sample dated polygon map series at x,y,t grid coordinates.    |
Raster     | Sample x,y raster at netcdf x,y grid coordinates              |


## **Create**

Create regular x,y,t netcdf file. 

Parameter | Description                                             |
--------- | ------------------------------------------------------  |
ncfile    | Output netcdf file path                                 |
raster    | Input raster file                                       |
start     | Start date YYYY-MM-DD                                   |
end       | End date YYYY-MM-DD                                     |
freq      | Time dimension frequency [DAILY,WEEKLY,MONTHLY,YEARLY]  |
tunit     | time units (default = "days since 1970-01-01 00:00:00") |
ttype     | time data type (default = "i4")                         |

## **Sample**

Sample dated polygon map series to netcdf grid using a query. 

Parameter  | Description                                     |
---------- | ----------------------------------------------- |
ncfile     | Netcdf file path                                |
varname    | Created variable name                           |
infeatures | Input polygon shapefile path                    |
target     | Extracted variable field name                   |
ineq       | Inequality used in query [<,<=,>=,>]            |
priority   | Priority field if multiple polygons match query |
datefield  | Input date field name                           |
dateformat | Input date format string (%Y-%m-%d)             |
dtype      | Created variable data type                      |
--revpri   | Reverse priority (-1 * priority)                |

Query 

    For each date in the time dimension
    select polygons where *polyDate* is *ineq* *gridDate*
    extract *targetfield* (choosing max *priorityfield* if n>1')

    
## **Raster**

Add a single raster variable to existing netcdf file

Parameter| Description           |
-------- | --------------------- |
ncfile   | Netcdf file path      | 
varname  | Created variable name |
raster   | Input raster          |
dtype    | Data type             |
