#!/usr/bin/python
#Copyright (C) 2013 Greg Miller <gmill002@gmail.com>
"""
netcdf

### Usage 
`ncgrid.py {Create,Sample,Raster} ...`

Subcommand | Description                                                           |
---------- | --------------------------------------------------------------------- |
Create     | Create regular x,y,t netcdf file from an input raster and date series |
Sample     | Sample dated polygon map series at x,y,t grid coordinates.            |
Raster     | Sample x,y raster at netcdf x,y grid coordinates                      |

"""

import arcpy
import argparse 
import datetime
from dateutil.rrule import rrule, YEARLY, MONTHLY, WEEKLY, DAILY
import netCDF4
import numpy as num
import os

arcpy.env.overwriteOutput = True

missvals = {'S1':'\x00','f4':9.96920996839e+36,'f8':9.96920996839e+36,'i1':-127,
            'i2':-32767,'i4':-2147483647,'i8': -9223372036854775806,'u1':255}

npytypes = {'S1':'?str?','f4':'float32','f8':'float64','i1':'int8','i2':'int16',
            'i4':'int32','i8':'int64','u1':'uint8'}

def Create(ncfile, raster, start, end, freq='YEARLY', ttype='i4', 
           tunit='days since 1970-01-01 00:00:00', xytype='i4'):
    """Create regular x,y,t netcdf file from an input raster and date series.
    
    Parameter  | Description                                      |
    ---------- | ------------------------------------------------ |
    ncfile     | output netcdf path                               |
    raster     | input raster path                                |
    start      | start date YYYY-MM-DD                            |
    end        | end date YYYY-MM-DD                              |
    freq       | sampling frequency [DAILY,WEEKLY,MONTHLY,YEARLY] | 

    """
    if not os.path.exists(os.path.join(os.path.dirname(ncfile), 'workspace')):
        os.makedirs(os.path.join(os.path.dirname(ncfile), 'workspace'))
    
    rootgrp = netCDF4.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
    desc = arcpy.Describe(raster)
    rootgrp.snapraster = raster
    rootgrp.spatialreference = desc.spatialReference.name
    rootgrp.cellsize = desc.meanCellHeight
    rootgrp.linear_unit = desc.spatialReference.linearUnitName

    xvals = (num.arange(desc.height) * desc.meanCellHeight) + desc.Extent.XMin + (0.5 * desc.meanCellHeight)    
    rootgrp.createDimension('x', len(xvals))    
    xcoord = rootgrp.createVariable('xcoord', xytype, ('x',))
    xcoord[:] = xvals.astype(npytypes[xytype])
    xcoord.units = 'x ' + desc.spatialReference.linearUnitName    

    yvals = (num.arange(desc.width) * desc.meanCellWidth) + desc.Extent.YMin + (0.5 * desc.meanCellWidth)    
    rootgrp.createDimension('y', len(yvals))
    ycoord = rootgrp.createVariable('ycoord', xytype, ('y',))
    ycoord[:] = yvals.astype(npytypes[xytype])
    ycoord.units = 'y ' + desc.spatialReference.linearUnitName    
    
    start, end = [datetime.datetime.strptime(t, "%Y-%m-%d") for t in [start, end]]
    dtvals = [dt for dt in rrule(eval(freq), dtstart=start, until=end)]
    timevals = [netCDF4.date2num(dt, tunit) for dt in dtvals]
    rootgrp.createDimension('t', len(timevals))
    tcoord = rootgrp.createVariable('tcoord', ttype, ('t',))
    tcoord[:] = num.array(timevals, dtype=npytypes[ttype])
    tcoord.units = tunit
    
    date = rootgrp.createVariable('date', 'S1', ('t',))
    date[:] = [datetime.datetime.strftime(dt, "%Y-%m-%d") for dt in dtvals]
        
    print rootgrp
    rootgrp.close()

def Raster(ncfile, varname, raster, dtype):
    """Add a single raster variable to existing netcdf file
    
    Parameter  | Description           |
    ---------- | --------------------- |
    ncfile     | netcdf file path      | 
    varname    | created variable name |
    raster     | input raster          |
    dtype      | data type             |
    
    """   
    rootgrp = netCDF4.Dataset(ncfile, 'a')
    arcpy.env.workspace = os.path.join(os.path.dirname(ncfile),'workspace')
    arcpy.CheckOutExtension("spatial")
    
    arcpy.ProjectRaster_management(raster, 'prj', rootgrp.snapraster, 'NEAREST', cell_size=rootgrp.cellsize)
    outExtractByMask = arcpy.sa.ExtractByMask('prj', rootgrp.snapraster)
    outExtractByMask.save('clip')
    
    var = rootgrp.createVariable(varname, dtype, ('x','y'), fill_value=missvals[dtype])
    var[:,:] = arcpy.RasterToNumPyArray('clip', nodata_to_value=missvals[dtype]).astype(npytypes[dtype])

    print var
    rootgrp.close()

def Sample(ncfile, varname, infeatures, target, ineq, priority, datefield, dateformat, dtype, revpri=False):
    """Sample dated polygon map series to netcdf grid using a query. 

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
    
    "revpri" should be used when one wants the earliest event
    
    """
    arcpy.env.snapraster = rootgrp.snapraster
    arcpy.env.workspace = os.path.join(os.path.dirname(ncfile),'workspace')
    
    # create new netcdf variable 
    rootgrp = netCDF4.Dataset(ncfile, 'a')
    var = rootgrp.createVariable(varname, dtype, ('x','y','t'), fill_value=missvals[dtype])
    
    # copy features to local file, convert date field 
    features = os.path.join(arcpy.env.workspace, os.path.basename(infeatures))
    arcpy.Copy_management(infeatures, features)
    ConvertDateField(features, datefield, dateformat)
    
    # Make new priority field with negative of original values
    # to select the earliest event when using days70 field
    if revpri: 
        if len(arcpy.ListFields(features, 'REVPRI')) == 0:
            arcpy.AddField_management(features, 'REVPRI', 'LONG')
        arcpy.CalculateField_management(features, 'REVPRI', '-{}'.format(priority))
        priority = 'REVPRI'
    
    arcpy.MakeFeatureLayer_management(features, "lyr")  
    for (i, date) in enumerate(rootgrp.variables['tcoord'][:]):
        # select all features meeting query criteria
        arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", '{} {} {}'.format('DAYS70', ineq, date))
        
        #workaround: sometimes the target raster is 'in use', so try another name        
        raster = os.path.join(arcpy.env.workspace, 'r{}'.format(i) )
        while 1: 
            try: 
                # extract to raster selecting by priority field
                arcpy.PolygonToRaster_conversion("lyr", target, raster, 'CELL_CENTER', priority, rootgrp.cellsize)
            except: 
                raster += 'x'
            else: 
                break
        
        var[:,:,i] = arcpy.RasterToNumPyArray(raster, nodata_to_value=missvals[dtype]).astype(npytypes[dtype])
    
    arcpy.Delete_management("lyr")
    print var
    rootgrp.close()
     
def ConvertDateField(inshape, infield, informat, newfield='DAYS70', newtype='LONG', outformat='days since 1970-01-01 00:00:00', dtype='i4'):   
    """Convert dates from an existing field/format into days since epoch integer field"""
    arcpy.AddField_management(inshape, newfield, newtype)    
    rows = arcpy.UpdateCursor(inshape)
    for row in rows:
        dt = datetime.datetime.strptime(str(getattr(row, infield)), informat)
        setattr(row, newfield, str(num.array(netCDF4.date2num(dt, outformat), dtype=npytypes[dtype])))
        rows.updateRow(row)

    
def Main():
    """command line parser with subcommands for create, raster, sample"""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="func")
    
    create = subparsers.add_parser('Create', help='Create a NetCDF grid from an input raster and date sequence. ')
    create.add_argument('ncfile', type=str, help='output netcdf path')
    create.add_argument('raster', type=str, help='input raster path')
    create.add_argument('start', type=str, help='start date YYYY-MM-DD')
    create.add_argument('end', type=str, help='end date YYYY-MM-DD')
    create.add_argument('freq', type=str, help='frequency of date samples',choices=['YEARLY','MONTHLY','WEEKLY','DAILY'])
    
    sample = subparsers.add_parser('Sample', help=("""Sample dated polygon map series to netcdf grid using a query. 
        For each date in the time dimension select polygons where <polyDate> is <ineq> <gridDate> extract <targetfield>
        choosing max <priorityfield> if n>1')"""))
    sample.add_argument('ncfile', type=str, help='netcdf file path')
    sample.add_argument('varname', type=str, help='created variable name')
    sample.add_argument('infeatures', type=str, help='input polygon shapefile path')
    sample.add_argument('target', type=str, help='extracted variable field name')
    sample.add_argument('ineq', type=str, help='inequality used in query')
    sample.add_argument('priority', type=str, help='priority field if multiple polygons meet query')
    sample.add_argument('datefield', type=str, help='input date field name')
    sample.add_argument('dateformat', type=str, help='input date format string (strptime format)')
    sample.add_argument('dtype', type=str, help='created variable data type')
    sample.add_argument('--revpri', default=False, action='store_true', help='reverse order indicated by priority field (-1 * x)')
    
    raster = subparsers.add_parser('Raster', help='Add a single raster variable to existing netcdf file')    
    raster.add_argument('ncfile', type=str, help='netcdf file path')
    raster.add_argument('varname', type=str, help='created variable name')
    raster.add_argument('raster', type=str, help='input raster')
    raster.add_argument('dtype', type=str, help='data type')    
    
    args = vars(parser.parse_args())
    FUNC = globals()[args.pop('func')]
    FUNC(**args)
    
if __name__ == '__main__':
    Main()
