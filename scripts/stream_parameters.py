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
    stream_segments = arcpy.da.UpdateCursor(hru.stream_path, ["FID", "TO_NODE", "tosegment", "k_coef", "x_coef"])
    compare_stream_segments = arcpy.da.SearchCursor(hru.stream_path, ["FID","FROM_NODE"])

    #Search all streams and find those whos from_nodes match another stream's to_node to determine tosegment param
    for  segment in stream_segments:
        to_node = segment[1] #to_node value
        #Check all other segments to see if they match the current segment's
        #This breaks as soon as a match is found therefore it is assumed that a stream does not split down stream
        for compare in compare_stream_segments:
            if to_node == compare[1]:  #compare to _node to from_node
                segment[2] = compare[0]+1 # tosegment = compare stream fid
                break 
            
        stream_segments.updateRow(segment)
        compare_stream_segments.reset()
     
    #Delete the structures created for generating the stream tosegment parameter
    del stream_segments, compare_stream_segments, compare, segment
 
    # Calculate the hru_segement by looking at how close points are to the streams. Take the smallest number as the answer
    logging.info("\nCalculating hru_segment")
    stream_segments = arcpy.da.SearchCursor(hru.stream_path, ["FID", "SHAPE@"])
    all_hrus = arcpy.da.UpdateCursor(hru.polygon_path, ["FID","HRU_SEG"])
    hru_centroids = arcpy.FeatureToPoint_management(hru.polygon_path)
    hru_centers = arcpy.da.SearchCursor(hru_centroids, ["FID", "SHAPE@X","SHAPE@Y"])

    hru_seg = {}

    for centroid in hru_centers:
        x=centroid[1]
        y=centroid[2]
        dist_dict = {}

        #Look at the combined distance a start and end point is away from the centroid
        #and use the min dist as the solution **Not perfect.        
        for strm in stream_segments:
            dist = []
            for pnt in [strm[1].firstPoint,strm[1].lastPoint]:
                dist.append((x-pnt.X)**2 + (y-pnt.Y)**2)

            dist_dict[strm[0]]=min(dist) 
  
        hru_seg[centroid[0]] = min(dist_dict,key = dist_dict.get)
        stream_segments.reset()
        
    #Assign the hru_segment.
    for hru in all_hrus:
        hru[1] = hru_seg[hru[0]]+1
        all_hrus.updateRow(hru)

    logging.info('\nDone!')
    del stream_segments,all_hrus,hru_centers,hru_seg,dist_dict

def average(lst):
    """
    Takes the average of a list
    """
    sum = 0
    for i in lst:
        sum+=i 
    return sum/len(lst)

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
