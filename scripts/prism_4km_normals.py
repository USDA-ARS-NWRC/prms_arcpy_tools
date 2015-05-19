#--------------------------------
# Name:         prism_4km_normals.py
# Purpose:      GSFLOW PRISM parameters from default 400m normals 
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2015-05-18
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import multiprocessing
import os
import re
import sys
##import tempfile
from time import clock

import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np

from support_functions import *

################################################################################

def gsflow_prism_parameters(config_path, data_name='ALL', 
                            overwrite_flag=False, debug_flag=False, ):
    """Calculate GSFLOW PRISM Parameters

    Args:
        config_file: Project config file path
        data_name -- the prism data type (ALL, PPT, TMAX, TMIN, etc.)
        ovewrite_flag: boolean, overwrite existing files
        debug_flag: boolean, enable debug level logging
    Returns:
        None
    """

    try:
        ## Initialize hru_parameters class
        hru = hru_parameters(config_path)

        ## Open input parameter config file
        inputs_cfg = ConfigParser.ConfigParser()
        try:
            inputs_cfg.readfp(open(config_path))
        except:
            logging.error('\nERROR: Config file could not be read, '+
                          'is not an input file, or does not exist\n'+
                          'ERROR: config_file = {0}\n').format(config_path)
            raise SystemExit()

        ## Log DEBUG to file
        log_file_name = 'gsflow_prism_normals_log.txt'
        log_console = logging.FileHandler(
            filename=os.path.join(hru.log_ws, log_file_name), mode='w')
        log_console.setLevel(logging.DEBUG)
        log_console.setFormatter(logging.Formatter('%(message)s'))
        logging.getLogger('').addHandler(log_console)
        logging.info('\nGSFLOW PRISM Parameters')

        ## PRISM
        prism_ws = inputs_cfg.get('INPUTS', 'prism_folder')
        prism_proj_method = inputs_cfg.get('INPUTS', 'prism_projection_method')
        prism_cs = inputs_cfg.getint('INPUTS', 'prism_cellsize')
        calc_prism_jh_coef_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_prism_jh_coef_flag')

        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(
                    hru.polygon_path))
            raise SystemExit()
        ## Check that PRISM folder is valid
        if not os.path.isdir(prism_ws):
            logging.error(
                '\nERROR: PRISM folder ({0}) does not exist'.format(prism_ws))
            raise SystemExit()
        proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if prism_proj_method.upper() not in proj_method_list:
            logging.error('\nERROR: PRISM projection method must be: {0}'.format(
                ', '.join(proj_method_list)))
            raise SystemExit()
        logging.debug('  Projection method:    {0}'.format(
            prism_proj_method.upper()))      

        ## Check other inputs
        if prism_cs <= 0:
            logging.error('\nERROR: PRISM cellsize must be greater than 0\n')
            raise SystemExit()

        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS 0'
        env.workspace = hru.param_ws
        env.scratchWorkspace = hru.scratch_ws      

        ## PRISM data names
        if data_name == 'ALL':
            data_name_list = ['PPT', 'TMAX', 'TMIN']
        else:
            data_name_list = [data_name]

        ## Set month list
        month_list = ['{0:02d}'.format(m) for m in range(1,13)]
        ##month_list.extend(['annual'])

        ## Check fields
        logging.info('\nAdding PRISM fields if necessary')
        for data_name in data_name_list:
            for month in month_list:
                add_field_func(
                    hru.polygon_path, '{0}_{1}'.format(data_name, month), 'DOUBLE')
            
        ## Process each PRISM data type
        logging.info('\nProjecting/clipping PRISM mean monthly rasters')
        for data_name in data_name_list:
            logging.info('\n{0}'.format(data_name))
            prism_normal_re = re.compile(
                'PRISM_(%s)_30yr_normal_4kmM2_(?P<month>\d{2})_bil.*$' % data_name,
                re.IGNORECASE)

            ## Search all files & subfolders in prism folder
            ##   for images that match data type
            input_raster_dict = dict()
            for root, dirs, files in os.walk(prism_ws):
                for file_name in files:
                    prism_normal_match = prism_normal_re.match(file_name)
                    if prism_normal_match and file_name.endswith('.bil'):
                        month_str = prism_normal_match.group('month')
                        input_raster_dict[month_str] = os.path.join(
                            prism_ws, root, file_name)
            if not input_raster_dict:
                logging.error(
                    ('\nERROR: No PRISM rasters were found matching the '+
                     'following pattern:\n  {0}\n\nDouble check that the script '+
                     'and folder are for the same resolution '+
                     '(800m vs 4km)\n\n').format(prism_normal_re.pattern))
                logging.error()
                raise SystemExit()

            ## PRISM input data workspace
            ##input_ws = os.path.join(prism_ws, data_name.lower())
            ##if not os.path.isdir(input_ws):
            ##    logging.error('\nERROR: The PRISM {0} folder does not exist'.format(
            ##        data_name.lower()))
            ##    raise SystemExit()

            ## PRISM output data workspace
            output_ws = os.path.join(
                hru.param_ws, data_name.lower()+'_rasters')
            if not os.path.isdir(output_ws):
                os.mkdir(output_ws)

            ## Remove all non year/month rasters in PRISM temp folder
            logging.info('  Removing existing PRISM files')
            for item in os.listdir(output_ws):
                if prism_normal_re.match(item):
                    os.remove(os.path.join(output_ws, item))

            ## Extract, project/resample, clip
            ## Process images by month
            zs_prism_dict = dict()
            ##env.extent = hru.extent
            for month in month_list:
                logging.info('  Month: {0}'.format(month))

                ## Projected/clipped PRISM raster
                input_raster = input_raster_dict[month]
                ##input_name = 'PRISM_{0}_30yr_normal_4kmM2_{1}_bil.bil'.format(
                ##    data_name.lower(), input_month)
                ##input_raster = os.path.join(input_ws, input_name)
                output_name = 'PRISM_{0}_30yr_normal_4kmM2_{1}.img'.format(
                    data_name.lower(), month)
                output_raster = os.path.join(output_ws, output_name)

                ## Set preferred transforms
                input_sr = Raster(input_raster).spatialReference
                transform_str = transform_func(hru.sr, input_sr)
                if transform_str:
                    logging.debug('  Transform: {0}'.format(transform_str))

                ## Project PRISM rasters to HRU coordinate system
                ## DEADBEEF - Arc10.2 ProjectRaster does not extent
                project_raster_func(
                    input_raster, output_raster, hru.sr,
                    prism_proj_method.upper(), prism_cs, transform_str,
                    '{0} {1}'.format(hru.ref_x, hru.ref_y), input_sr, hru)
                ##arcpy.ProjectRaster_management(
                ##    input_raster, output_raster, hru.sr,
                ##    prism_proj_method.upper(), prism_cs, transform_str,
                ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
                ##    input_sr)

                ## Save parameters for calculating zonal stats
                zs_field = '{0}_{1}'.format(data_name, month)
                zs_prism_dict[zs_field] = [output_raster, 'MEAN']

                ## Cleanup
                del input_raster, output_raster, output_name
                del input_sr, transform_str, zs_field

            ## Cleanup
            ##arcpy.ClearEnvironment('extent')

            ## Calculate zonal statistics
            logging.info('\nCalculating PRISM zonal statistics')
            zonal_stats_func(zs_prism_dict, hru.polygon_path, hru.point_path, hru)
            del zs_prism_dict
        
        ## Jensen-Haise Potential ET air temperature coefficient
        ## Update Jensen-Haise PET estimate using PRISM air temperature
        ## DEADBEEF - First need to figure out month with highest Tmax
        ##            Then get Tmin for same month
        if calc_prism_jh_coef_flag:
            logging.info('\nRe-Calculating JH_COEF_HRU')
            logging.info('  Using PRISM temperature values')
            tmax_field_list = ['!TMAX_{0:02d}!'.format(m) for m in range(1,13)]
            tmin_field_list = ['!TMIN_{0:02d}!'.format(m) for m in range(1,13)]           
            tmax_expr = 'max([{0}])'.format(','.join(tmax_field_list))
            arcpy.CalculateField_management(
                hru.polygon_path, hru.jh_tmax_field, tmax_expr, 'PYTHON')
            ## Sort TMAX and get TMIN for same month
            tmin_expr = 'max(zip([{0}],[{1}]))[1]'.format(
                ','.join(tmax_field_list), ','.join(tmin_field_list))
            arcpy.CalculateField_management(
                hru.polygon_path, hru.jh_tmin_field, tmin_expr, 'PYTHON')    

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

    finally:
        try: arcpy.CheckInExtension('Spatial')
        except: pass
        ##arcpy.ResetEnvironments()

################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='PRISM 4km Normals',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
    parser.add_argument(
        '--type', default='ALL', 
        help='PRISM Data Type (TMAX, TMIN, PPT, ALL)')
    parser.add_argument(
        '-o', '--overwrite', default=False, action="store_true", 
        help='Force overwrite of existing files')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()

    ## Create Basic Logger
    logging.basicConfig(level=args.loglevel, format='%(message)s')

    #### Get GSFLOW config file
    ##ini_re = re.compile('\w*.ini$', re.I)
    ##try: 
    ##    ini_path = sys.argv[1]
    ##except IndexError:
    ##    ini_path = get_ini_file(workspace, ini_re, 'gsflow_prism_parameters')
    ##del ini_re

    ## Get PRISM data name as argument only if a config file was also set
    ##try: 
    ##    data_name = sys.argv[2]
    ##except IndexError:
    ##    data_name = get_prism_data_name()
    ##if data_name not in ['PPT', 'TMAX', 'TMIN']:
    ##    data_name = 'ALL'

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', os.getcwd()))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    ## Calculate GSFLOW PRISM Parameters
    gsflow_prism_parameters(
        config_path=args.ini, data_name=args.type, 
        overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
