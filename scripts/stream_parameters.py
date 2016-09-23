#--------------------------------
# Name:         stream_parameters.py
# Purpose:      GSFLOW stream parameters
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2016-02-26
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import os
# import re
import shutil
import subprocess
import sys
# from time import clock, sleep

import arcpy
from arcpy import env
from arcpy.sa import *

# import numpy as np

from support_functions import *
from matplotlib.lines import segment_hits


def stream_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate PRMS Stream Parameters

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    # Initialize hru_parameters class
    hru = HRUParameters(config_path)

    # Log DEBUG to file
    log_file_name = 'stream_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nPRMS Stream Parameters')


    # check the polygon path
    hru.check_polygon_path()
    
    # read the stream parameters
    hru.read_stream_parameters()
    
    # Input folders
    stream_temp_ws = os.path.join(hru.param_ws, 'stream_rasters')
    if not os.path.isdir(stream_temp_ws):
        os.mkdir(stream_temp_ws)
          
    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    # env.pyramid = 'PYRAMIDS -1'
    env.pyramid = 'PYRAMIDS 0'
    env.workspace = stream_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Add fields if necessary to the HRU
    logging.info('\nAdding fields if necessary')
    add_field_func(hru.polygon_path, hru.hru_segment, 'LONG')
    add_field_func(hru.stream_path, hru.k_coef, 'DOUBLE')
    add_field_func(hru.stream_path, hru.obsin_segment, 'LONG')
    add_field_func(hru.stream_path, hru.tosegment, 'LONG')
    add_field_func(hru.stream_path, hru.x_coef, 'DOUBLE')
    
    # Calculate the TOSEGMENT, k_coef, x_coef
    logging.info("\nCalculating tosegment, k_coef, and x_coef parameters")
    stream_segments = arcpy.da.UpdateCursor(hru.stream_path, ["OBJECTID", "to_node", "tosegment", "k_coef", "x_coef"])
    compare_stream_segments = arcpy.da.SearchCursor(hru.stream_path, ["OBJECTID","from_node"])

    #Search all streams and find those whos from_nodes match another stream's to_node to determine tosegment param
    for  segment in stream_segments:
        to_node = segment[1] #to_node value
        #Check all other segments to see if they match the current segment's
        #This breaks as soon as a match is found therefore it is assumed that a stream does not split down stream
        for compare in compare_stream_segments:
            if to_node == compare[1]:  #compare to _node to from_node
                segment[2] = compare[0] # tosegment = compare stream objectid
                break 
            
        stream_segments.updateRow(segment)
        compare_stream_segments.reset()
     
    #Delete the structures created for generating the stream tosegment parameter
    del stream_segments, compare_stream_segments, compare, segment
 
    # Calculate the hru_segement
    logging.info("\nCalculating hru_segment")
    stream_segments = arcpy.da.SearchCursor(hru.stream_path, ["OBJECTID", "grid_code"])
    all_hrus = arcpy.da.UpdateCursor(hru.polygon_path, ["OBJECTID","grid_code","HRU_SEG"])

    #Search all streams and HRU to find matching Gridcode to assign the HRU_segment param
    for  hru in all_hrus:
        hru_gc = hru[1] #to_node value
        #Check all other segments to see if they match the current segment's
        #This breaks as soon as a match is found therefore it is assumed that a stream does not split down stream
        for stream in stream_segments:
            if hru_gc == stream[1]:  #compare to _node to from_node
                hru[2] = stream[0]
                break 
        all_hrus.updateRow(hru)
        stream_segments.reset()
    
    # Get stream length for each cell
#     logging.info("Stream length")
#     arcpy.MakeFeatureLayer_management(hru.polygon_path, hru_polygon_lyr)
#     arcpy.SelectLayerByAttribute_management(
#         hru_polygon_lyr, "NEW_SELECTION",
#         ' \"{0}\" = 1 And "{1}" != 0'.format(hru.type_field, hru.iseg_field))
#     length_path = os.path.join('in_memory', 'length')
#     arcpy.Intersect_analysis(
#         [hru_polygon_lyr, hru.streams_path],
#         length_path, "ALL", "", "LINE")
#     arcpy.Delete_management(hru_polygon_lyr)
#      
#     
#     arcpy.CalculateField_management(
#         hru.stream_path, length_field, '!shape.length@meters!', "PYTHON")
#     length_dict = defaultdict(int)
#      
#         # DEADBEEF - This probably needs a maximum limit
#     for row in arcpy.da.SearchCursor(
#         length_path, [hru.id_field, length_field]):
#         length_dict[int(row[0])] += int(row[1])
#     fields = [hru.type_field, hru.iseg_field, hru.rchlen_field, hru.id_field]
#     with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
#         for row in update_c:
#             if (int(row[0]) == 1 and int(row[1]) != 0):
#                 row[2] = length_dict[int(row[3])]
#             else:
#                 row[2] = 0
#             update_c.updateRow(row)
#     del length_dict, length_field, fields, hru_polygon_lyr

#     # Set environment parameters
#     env.extent = hru.extent
#     env.cellsize = hru.cs
#     env.outputCoordinateSystem = hru.sr
# 
#     # Build rasters
#     if output_rasters_flag:
#         logging.info("\nOutput model grid rasters")
#         arcpy.PolygonToRaster_conversion(
#             hru.polygon_path, hru.type_field, hru_type_raster,
#             "CELL_CENTER", "", hru.cs)
#         arcpy.PolygonToRaster_conversion(
#             hru.polygon_path, hru.dem_adj_field, dem_adj_raster,
#             "CELL_CENTER", "", hru.cs)
#         arcpy.PolygonToRaster_conversion(
#             hru.polygon_path, hru.iseg_field, iseg_raster,
#             "CELL_CENTER", "", hru.cs)
#         arcpy.PolygonToRaster_conversion(
#             hru.polygon_path, hru.irunbound_field, irunbound_raster,
#             "CELL_CENTER", "", hru.cs)
#         arcpy.PolygonToRaster_conversion(
#             hru.polygon_path, hru.segbasin_field, segbasin_raster,
#             "CELL_CENTER", "", hru.cs)
#         arcpy.PolygonToRaster_conversion(
#             hru.polygon_path, hru.subbasin_field, subbasin_raster,
#             "CELL_CENTER", "", hru.cs)
# 
#     # Build rasters
#     if output_ascii_flag:
#         logging.info("Output model grid ascii")
#         arcpy.RasterToASCII_conversion(hru_type_raster, hru_type_ascii)
#         arcpy.RasterToASCII_conversion(dem_adj_raster, dem_adj_ascii)
#         arcpy.RasterToASCII_conversion(iseg_raster, iseg_ascii)
#         arcpy.RasterToASCII_conversion(irunbound_raster, irunbound_ascii)
#         arcpy.RasterToASCII_conversion(segbasin_raster, segbasin_ascii)
#         arcpy.RasterToASCII_conversion(subbasin_raster, subbasin_ascii)
#         sleep(5)
    logging.info('\nDone!')

def cell_distance(cell_a, cell_b, cs):
    """"""
    ai, aj = cell_a
    bi, bj = cell_b
    return math.sqrt((ai - bi) ** 2 + (aj - bj) ** 2) * cs

# def calc_stream_width(flow_acc):
#    return -2E-6 * flow_acc ** 2 + 0.0092 * flow_acc + 1


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Stream Parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
    parser.add_argument(
        '-o', '--overwrite', default=False, action="store_true",
        help='Force overwrite of existing files')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()

    # Convert input file to an absolute path
    if os.path.isfile(os.path.abspath(args.ini)):
        args.ini = os.path.abspath(args.ini)
    return args


if __name__ == '__main__':
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format('Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', os.getcwd()))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    # Calculate PRMS Muskingum Stream Parameters
    stream_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
