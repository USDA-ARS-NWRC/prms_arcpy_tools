#--------------------------------
# Name:        run_all
# Purpose:     PRMS Execute All
# Notes:        ArcGIS 10.2 Version
# Author:       Micah Johnson
# Created       2016-09-13
# Python:       2.7
#--------------------------------

import argparse
import ConfigParser
import datetime as dt
import logging
import os
import sys

import arcpy
from arcpy import env
from arcpy.sa import *

from support_functions import *
from hru_parameters import hru_parameters
from dem_parameters import dem_parameters
from veg_parameters import veg_parameters
from soil_raster_prep import soil_raster_prep
from soil_parameters import soil_parameters
from impervious_parameters import impervious_parameters
from prism_4km_normals import prism_4km_parameters
from ppt_ratio_parameters import ppt_ratio_parameters
from stream_parameters  import stream_parameters
from prms_template_fill import prms_template_fill

def calculate_all_parameters(config_path, data_name='ALL', overwrite_flag=False, debug_flag=False):
    """
    Calculate all PRMS Parameters

    Executes all the parameter scripts in order
    to build a parameter file for PRMS

    The execution order is as follows:
    hru_parameters
    dem_parameters
    veg_parameters
    soil_raster_prep
    soil_parameters
    impervious_parameters
    prism_4km_normals  / prism_800m_normals
    ppt_ratio_parameters
    stream_parameters
    prms_template_fill

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """
    logging.info('\nCalculating PRMS HRU Parameters...')
    hru_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)

    logging.info("\nFinished!")

    logging.info("\nCalculating PRMS DEM Parameters...")
    dem_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\nCalculating PRMS Vegetation Parameters...")
    veg_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\nPreparing Soil Rasters...")
    soil_raster_prep(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\nCalculating PRMS Soil Parameters...")
    soil_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\nCalculating PRMS Impervious Parameters...")
    impervious_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\nCalculating PRISM 4Km Parameters...")
    prism_4km_parameters(
        config_path=args.ini, data_name=args.type,
        overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\nCalculating PPT Ratio Parameters...")
    ppt_ratio_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\nCalculating PRMS Stream Parameters...")
    stream_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nFinished!")

    logging.info("\n Writing Parameters to Input File for PRMS...")
    prms_template_fill(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
    logging.info("\nParameters are now written to file and can be used for PRMS Simulations.")
def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='PRMS Template Fill',
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
    parser.add_argument(
        '--type', default='ALL',
        help='PRISM Data Type (TMAX, TMIN, PPT, ALL)')

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

    # Calculate PRMS Parameters
    calculate_all_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
