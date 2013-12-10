ncgrid
======

Interpolate dated polygon map series to a 3-d (x,y,t) lattice (netcdf file)


Usage
-----

`ncgrid.py {Create,Sample,Raster} [command options]`

==========  ==============================================================                                                                            
Command     Description                                                   
==========  ==============================================================                                                                            
create      Create x,y,t netcdf file from an input raster and date series 
sample      Sample dated polygon map series at x,y,t grid coordinates.    
raster      Sample x,y raster at netcdf x,y grid coordinates              
==========  ==============================================================                                                                            


Create
------

Create x,y,t lattice netcdf file. 

=========  ==============================================================                                                                            
Parameter  Description                                             
=========  ==============================================================                                                                            
ncfile     Created netcdf file path                                 
raster     Input raster file. Defines xy lattice                                      
start      Start date YYYY-MM-DD                                   
end        End date YYYY-MM-DD                                     
freq       Time dimension frequency [DAILY,WEEKLY,MONTHLY,YEARLY]  
tunit      time units (default = "days since 1900-01-01 00:00:00") 
ttype      time data type (default = "i4")                         
=========  ==============================================================                                                                            

sample
------

sample dated polygon map series to netcdf grid using a query. 

==========  ==============================================================                                                                            
Parameter   Description                                     
==========  ==============================================================                                                                            
ncfile      Netcdf file path                                
varname     Created variable name                           
infeatures  Input polygon shapefile path                    
target      Extracted variable field name                   
ineq        Inequality used in query [<,<=,>=,>]            
priority    Priority field if multiple polygons match query 
datefield   Input date field name                           
dateformat  Input date format string (%Y-%m-%d)             
dtype       Created variable data type                      
--revpri    Reverse priority (-1 * priority)                
==========  ==============================================================                                                                            

Query
::
    
    For each date in the time dimension
    select polygons where *polyDate* is *ineq* *gridDate*
    extract *targetfield* (choosing max *priorityfield* if n>1')

option `--revpri` should be used in queries of the 'next' feature of interest 
(earliest feature that is greater than current date). 

    
raster
------

Add a single raster variable to existing netcdf file

=========  ==============================================================                                                                            
Parameter  Description           
=========  ==============================================================                                                                            
ncfile     Netcdf file path       
varname    Created variable name 
raster     Input raster          
dtype      Data type             
=========  ==============================================================                                                                            
