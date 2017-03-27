"""
The parameters sometimes to not have default values for larger areas. This
is typically seen from the NRCS soils database. So take the polygons that do
not have values and estimate based on their neighbors.

20170327 Scott Havens
"""

import argparse
import datetime as dt
import logging
import os
import sys
# from time import clock, sleep

import arcpy
from arcpy import env

import numpy as np

def count_selected(layer):
    result = arcpy.GetCount_management(layer)
    selected = int(result.getOutput(0))
    logging.debug('Selected {} features'.format(selected))
    
    return selected

def estimate(file_name, field, null_value):
    """
    Estimate the parameters, 
    """
    
    # open the shapefile
    f = arcpy.da.UpdateCursor(file_name, ["FID","SHAPE@",field],
                              where_clause='{}={}'.format(field, null_value))
        
    # Make a layer from the feature class
    arcpy.MakeFeatureLayer_management(file_name, "lyr") 
    ret = True
    
    for h in f:
        
        logging.debug('Feature FID = {}'.format(h[0]))
                        
        # create the initial selection
        arcpy.SelectLayerByAttribute_management("lyr", 'NEW_SELECTION', 'FID={}'.format(h[0]))
        count_selected("lyr")
        
        # select where the boundary touches
        arcpy.SelectLayerByLocation_management("lyr", 'BOUNDARY_TOUCHES', "lyr")
        selected = count_selected("lyr")
                
        # get the selected values for sand
        fld = np.zeros((selected,))
        area = np.zeros((selected,))
        for i,row in enumerate(arcpy.SearchCursor("lyr")):
            fld[i] = row.getValue(field)
            area[i] = row.Shape.area
            
        # remove any null values
        idx = fld == null_value
        fld = fld[~idx]
        area = area[~idx]

        # update with the average
        w = area/np.sum(area)
        
        if np.sum(w) == 0:
            ret = False
        else:
            h[2] = np.average(fld, weights=w)
            f.updateRow(h)
        
            logging.debug('Updated value to {}'.format(h[2]))
        
    return ret
    

def parameter_estimator(file_name, field, null_value=0):
    """
    Estimate the parameters for the given shapfile
    """
    
    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    # env.pyramid = 'PYRAMIDS -1'
    env.pyramid = 'PYRAMIDS 0'
#     env.workspace = stream_temp_ws
    env.scratchWorkspace = 'in_memory'
    
    for it in range(3):
        logging.debug('Iteration {}'.format(it))
        
        ret = estimate(file_name, field, null_value)
        
        if ret:
            break

    logging.info('\nDone!')


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Parameter Estimator!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-f', '--file', required=True,
        help='Input shapefile')
    parser.add_argument(
        '-c', '--field', required=True,
        help='Shapefile field name')
    parser.add_argument(
        '-n', '--null',
        help='Null value to overwrite, default 0')
    args = parser.parse_args()

    # Convert input file to an absolute path
    if os.path.isfile(os.path.abspath(args.file)):
        args.file = os.path.abspath(args.file)
    if args.null is None:
        args.null = 0
    else:
        args.null = int(args.null)
    return args


if __name__ == '__main__':
    args = arg_parse()

    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format('Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', os.getcwd()))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    # Calculate PRMS Muskingum Stream Parameters
    parameter_estimator(
        file_name=args.file,
        field=args.field,
        null_value=args.null)
    
    