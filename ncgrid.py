"""
ncgrid
======

Interpolate dated polygon map series to a 3-d (x,y,t) lattice (netcdf file)
"""
import argparse 
import datetime
from dateutil.rrule import rrule, YEARLY, MONTHLY, WEEKLY, DAILY
import os

import arcpy
import netCDF4
import numpy as num
from pathlib import Path

arcpy.env.overwriteOutput = True

missvals = {'S1':'\x00','f4':9.96920996839e+36,'f8':9.96920996839e+36,'i1':-127,
            'i2':-32767,'i4':-2147483647,'i8': -9223372036854775806,'u1':255}

npytypes = {'S1':'?str?','f4':'float32','f8':'float64','i1':'int8','i2':'int16',
            'i4':'int32','i8':'int64','u1':'uint8'}
    
def create(ncfile, raster, start, end, freq='YEARLY', ttype='i4', 
           tunit='days since 1970-01-01 00:00:00', xytype='i4'):
    """Create a NetCDF file with dimensions based on input raster and frequency.
    
    Parameter  Description  
    =========  ===================================
    ncfile     created netcdf file path
    raster     ESRI GRID raster. cell centers define the lattice points. linear_unit in meters. 
    start      timeseries start date "YYYY-MM-DD"             
    end        timeseries end date "YYYY-MM-DD"
    freq       sampling frequency [YEARLY, MONTHLY, WEEKLY, DAILY]
    ttype      time datatype
    tunit      time units
    xytype     xy units

    """
    ws = Path(ncfile).parent['workspace']
    if not ws.exists():
        ws.mkdir()
    
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

def raster(ncfile, varname, raster, dtype):
    """Add a raster xy variable to existing netcdf file. 
    
    Input raster is automatically projected and resampled at the input netcdf xy coordinates. 
    
    """   
    rootgrp = netCDF4.Dataset(ncfile, 'a')
    arcpy.env.workspace = str(Path(ncfile).parent()['workspace'])
    #arcpy.env.workspace = os.path.join(os.path.dirname(ncfile),'workspace')
    arcpy.CheckOutExtension("spatial")
    
    arcpy.ProjectRaster_management(raster, 'prj', rootgrp.snapraster, 'NEAREST', cell_size=rootgrp.cellsize)
    outExtractByMask = arcpy.sa.ExtractByMask('prj', rootgrp.snapraster)
    outExtractByMask.save('clip')
    
    var = rootgrp.createVariable(varname, dtype, ('x','y'), fill_value=missvals[dtype])
    var[:,:] = arcpy.RasterToNumPyArray('clip', nodata_to_value=missvals[dtype]).astype(npytypes[dtype])

    print var
    rootgrp.close()

def sample(ncfile, varname, infeatures, target, ineq, priority, datefield, dateformat, dtype, revpri=False):
    """Sample a dated polygon map series on the input NetCDF grid. 
    
    For each date coordinate in the netcdf time dimension 
    
    select polygons where <polyDate> is <ineq> <gridDate> extract <targetfield> 
    choosing max <priorityfield> if n>1
    
    If interpolation is last ( < | <= ) then priority field DAYS70 is fine. 
    If interpolation is next ( > | >= ) then revpri must be TRUE or input an inverted numeric date field
    
    """
    # open the dataset, set env variables 
    rootgrp = netCDF4.Dataset(ncfile, 'a')
    arcpy.env.snapraster = rootgrp.snapraster
    #arcpy.env.workspace = os.path.join(os.path.dirname(ncfile),'workspace')
    arcpy.env.workspace = str(Path(ncfile).parent()['workspace'])
    
    # copy input shapefile and convert dates 
    temp = os.path.join(arcpy.env.workspace, os.path.basename(infeatures))
    arcpy.Copy_management(infeatures, temp)
    convert_shapefile_date_field(temp, datefield, dateformat)
    
    # reverse priority: create new field that is the reverse of the priority field
    if revpri:
        if len(arcpy.ListFields(temp, 'REVPRI')) == 0:
            arcpy.AddField_management(temp, 'REVPRI', 'LONG')
    
        arcpy.CalculateField_management(temp, 'REVPRI', '-1 * !{}!'.format(priority), 'PYTHON')
        priority = 'REVPRI'
    
    # create in memory polygons layer
    arcpy.MakeFeatureLayer_management(temp, "lyr")  
    
    #create output netcdf variable
    var = rootgrp.createVariable(varname, dtype, ('x','y','t'), fill_value=missvals[dtype])
    
    # for each date: select polygons and extract to a raster then to the output variable slice  
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
     
def convert_shapefile_date_field(inshape, infield, informat, newfield='DAYS70', newtype='LONG', outformat='days since 1970-01-01 00:00:00', dtype='i4'):   
    """Convert dates from an existing field/format into days since epoch integer field"""
    arcpy.AddField_management(inshape, newfield, newtype)    
    rows = arcpy.UpdateCursor(inshape)
    for row in rows:
        dt = datetime.datetime.strptime(str(getattr(row, infield)), informat)
        setattr(row, newfield, str(num.array(netCDF4.date2num(dt, outformat), dtype=npytypes[dtype])))
        rows.updateRow(row)

    
def command_args():
    """command line parser with subcommands for create, raster, sample"""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    create = subparsers.add_parser('create', help='Create a NetCDF grid from an input raster and date sequence. ')
    create.add_argument('ncfile', type=str, help='output netcdf path')
    create.add_argument('raster', type=str, help='input raster path')
    create.add_argument('start', type=str, help='start date YYYY-MM-DD')
    create.add_argument('end', type=str, help='end date YYYY-MM-DD')
    create.add_argument('freq', type=str, help='frequency of date samples',choices=['YEARLY','MONTHLY','WEEKLY','DAILY'])

    sample = subparsers.add_parser('sample', help=("""Sample dated polygon map series to netcdf grid using a query. 
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

    raster = subparsers.add_parser('raster', help='Add a single raster variable to existing netcdf file')
    raster.add_argument('ncfile', type=str, help='netcdf file path')
    raster.add_argument('varname', type=str, help='created variable name')
    raster.add_argument('raster', type=str, help='input raster')
    raster.add_argument('dtype', type=str, help='data type')
    
    return vars(parser.parse_args())
    
if __name__ == '__main__':
    args = command_args()
    globals()[args.pop('cmd')](**args)
