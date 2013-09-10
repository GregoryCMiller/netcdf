"""
netcdf
======

Interpolate dated polygon map series to a regular 3-d grid (x,y,t) stored as a netcdf variable

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
    """Create a NetCDF grid from an input raster and date sequence parameters.
    
    Create a NETCDF3 CLASSIC file that has dimensions (x,y,t). 
    
    
    Dimensions x,y 
        - Coordinates are defined as input raster cell centers. 
        - linear unit must be meters. decimal degrees might cause problems
    
    
    Dimension t is composed of <count> equally spaced integer dates from <start> to <end> inclusive. 
    
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
    """Add a raster xy variable to existing netcdf file. 
    
    Input raster is automatically projected and resampled at the input netcdf xy coordinates. 
    
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
    """Sample a dated polygon map series on the input NetCDF grid. 
    
    For each date coordinate in the netcdf time dimension 
    
    select polygons where <polyDate> is <ineq> <gridDate> extract <targetfield> 
    choosing max <priorityfield> if n>1
    
    If interpolation is last ( < | <= ) then priority field DAYS70 is fine. 
    If interpolation is next ( > | >= ) then revpri must be TRUE or input an inverted numeric date field
    
    """
    
    rootgrp = netCDF4.Dataset(ncfile, 'a')
    arcpy.env.snapraster = rootgrp.snapraster
    arcpy.env.workspace = os.path.join(os.path.dirname(ncfile),'workspace')
    
    temp = os.path.join(arcpy.env.workspace, os.path.basename(infeatures))
    arcpy.Copy_management(infeatures, temp)
    ConvertDateField(temp, datefield, dateformat)
    
    if revpri: # option to make new priority field with negative of original values
        if len(arcpy.ListFields(temp, 'REVPRI')) == 0:
            arcpy.AddField_management(temp, 'REVPRI', 'LONG')
    
        arcpy.CalculateField_management(temp, 'REVPRI', '-{}'.format(priority))
        priority = 'REVPRI'
        
    arcpy.MakeFeatureLayer_management(temp, "lyr")  
    
    var = rootgrp.createVariable(varname, dtype, ('x','y','t'), fill_value=missvals[dtype])
    
    for (i, date) in enumerate(rootgrp.variables['tcoord'][:]):
        arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", '{} {} {}'.format('DAYS70', ineq, date))
        raster = os.path.join(arcpy.env.workspace, 'r{}'.format(i) )
        while 1: #workaround: sometimes the target raster is 'in use', so try another name
            try: 
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
    globals()[args.pop('func')](**args)
    
if __name__ == '__main__':
    Main()
