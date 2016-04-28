#--------------------------------
# Name:         impervious_parameters.py
# Purpose:      GSFLOW impervious parameters
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2016-02-26
# Python:       2.7
#--------------------------------

import argparse
# from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import os
# import re
import sys

import arcpy
from arcpy import env
from arcpy.sa import *
# import numpy as np

from support_functions import *


def impervious_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate PRMS Impervious Parameters

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
    log_file_name = 'impervious_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nPRMS Impervious Parameters')

    # Check input paths
    hru.check_polygon_path()

    # Check impervious parameters from config file
    hru.read_impervious_parameters()     
   

    # Build output folder if necessary
    imperv_temp_ws = os.path.join(hru.param_ws, 'impervious_rasters')
    if not os.path.isdir(imperv_temp_ws):
        os.mkdir(imperv_temp_ws)
    # Output paths
    imperv_path = os.path.join(imperv_temp_ws, 'impervious_cover.img')


    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    env.workspace = imperv_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Check field
    logging.info('\nAdding impervious fields if necessary')
    add_field_func(hru.polygon_path, hru.imperv_pct_field, 'DOUBLE')
    # add_field_func(hru.polygon_path, hru.carea_min_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.carea_max_field, 'DOUBLE')

    # Available Water Capacity (AWC)
    logging.info('\nProjecting/clipping impervious cover raster')
    imperv_orig_sr = Raster(hru.imperv_orig_path).spatialReference
    logging.debug('  Impervious GCS:  {0}'.format(
        imperv_orig_sr.GCS.name))
    # Remove existing projected raster
    if arcpy.Exists(imperv_path):
        arcpy.Delete_management(imperv_path)
    # Set preferred transforms
    transform_str = transform_func(hru.sr, imperv_orig_sr)
    logging.debug('  Transform: {0}'.format(transform_str))
    logging.debug('  Projection method: NEAREST')
    # Project impervious raster
    # DEADBEEF - Arc10.2 ProjectRaster does not extent
    # env.extent = hru.extent
    project_raster_func(
        hru.imperv_orig_path, imperv_path, hru.sr,
        hru.imperv_proj_method, hru.imperv_cs, transform_str,
        imperv_orig_sr, hru)
    # arcpy.ProjectRaster_management(
    #    imperv_orig_path, imperv_path, hru.sr,
    #    imperv_proj_method, imperv_cs, transform_str,
    #    '{0} {1}'.format(hru.ref_x, hru.ref_y),
    #    imperv_orig_sr)
    # arcpy.ClearEnvironment('extent')

    # List of rasters, fields, and stats for zonal statistics
    zs_imperv_dict = dict()
    zs_imperv_dict[hru.imperv_pct_field] = [imperv_path, 'MEAN']
    # zs_imperv_dict[hru.carea_min_field] = [imperv_path, 'MEAN']
    # zs_imperv_dict[hru.carea_max_field] = [imperv_path, 'MEAN']

    # Calculate zonal statistics
    logging.info('\nCalculating zonal statistics')
    zonal_stats_func(zs_imperv_dict, hru.polygon_path, hru)

    # Calculate CAREA_MIN / CAREA_MAX
    logging.info('\nCalculating CAREA_MIN / CAREA_MAX')
    if hru.imperv_pct_flag:
        arcpy.CalculateField_management(
            hru.polygon_path, hru.imperv_pct_field,
            '0.01 * !{0}!'.format(hru.imperv_pct_field), 'PYTHON')
        # arcpy.CalculateField_management(
        #    hru.polygon_path, hru.carea_min_field,
        #    '0.01 * !{0}!'.format(hru.imperv_pct_field), 'PYTHON')
        arcpy.CalculateField_management(
            hru.polygon_path, hru.carea_max_field,
            '0.01 * !{0}!'.format(hru.imperv_pct_field), 'PYTHON')
    else:
        # arcpy.CalculateField_management(
        #    hru.polygon_path, hru.carea_min_field,
        #    '!{0}!'.format(hru.imperv_pct_field), 'PYTHON')
        arcpy.CalculateField_management(
            hru.polygon_path, hru.carea_max_field,
            '!{0}!'.format(hru.imperv_pct_field), 'PYTHON')
        
    logging.info('Done!')


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Impervious Cover Parameters',
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

    # Calculate GSFLOW Impervious Parameters
    impervious_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
