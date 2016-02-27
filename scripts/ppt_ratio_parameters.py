#--------------------------------
# Name:         ppt_ratio_parameters.py
# Purpose:      GSFLOW PPT ratio parameters
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
# import multiprocessing
import os
# import re
import sys
# from time import clock

import arcpy
from arcpy import env
from arcpy.sa import *
# import numpy as np

import support_functions


def ppt_ratio_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate GSFLOW PPT Ratio Parameters

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    # Initialize hru_parameters class
    hru = support_functions.HRUParameters(config_path)

    # Open input parameter config file
    inputs_cfg = ConfigParser.ConfigParser()
    try:
        inputs_cfg.readfp(open(config_path))
    except:
        logging.error('\nERROR: Config file could not be read, ' +
                      'is not an input file, or does not exist\n' +
                      'ERROR: config_file = {0}\n').format(config_path)
        sys.exit()

    # Log DEBUG to file
    log_file_name = 'ppt_ratio_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nGSFLOW PPT Ratio Parameters')

    # Units
    ppt_obs_units = support_functions.get_param(
        'ppt_obs_units', 'mm', inputs_cfg).lower()
    ppt_units_list = ['mm', 'cm', 'm', 'in', 'ft']
    # Compare against the upper case of the values in the list
    #   but don't modify the acceptable units list
    if ppt_obs_units not in ppt_units_list:
        logging.warning(
            ('WARNING: Invalid PPT obs. units ({0})\n  ' +
             'Valid units are: {1}').format(
                ppt_obs_units, ', '.join(ppt_units_list)))

    # PPT Zones
    set_ppt_zones_flag = inputs_cfg.getboolean('INPUTS', 'set_ppt_zones_flag')
    if set_ppt_zones_flag:
        ppt_zone_orig_path = inputs_cfg.get('INPUTS', 'ppt_zone_path')
        ppt_zone_field = inputs_cfg.get('INPUTS', 'ppt_zone_field')
    # If a zone shapefile is not used, PPT must be set manually
    else:
        ppt_obs_list = inputs_cfg.get('INPUTS', 'ppt_obs_list')
        try:
            ppt_hru_id = inputs_cfg.getint('INPUTS', 'ppt_hru_id')
        except:
            ppt_hru_id = 0
        # Check that values are floats
        try:
            ppt_obs_list = map(float, ppt_obs_list.split(','))
        except ValueError:
            logging.error(
                '\nERROR: ppt_obs_list (mean monthly precipitation) ' +
                'values could not be parsed as floats')
            sys.exit()

        # Check that there are 12 values
        if len(ppt_obs_list) != 12:
            logging.error(
                '\nERROR: There must be exactly 12 mean monthly ' +
                'observed precipitation values based to ppt_obs_list')
            sys.exit()

        # Check that HRU_ID is valid
        logging.info('  PPT HRU_ID: {0}'.format(ppt_hru_id))
        arcpy.MakeTableView_management(
            hru.polygon_path, "test_layer",
            "{0} = {1}".format(hru.id_field, ppt_hru_id))
        if ppt_hru_id == 0:
            logging.info(
                '    Assuming ppt_ratios should not be forced to 1 at a cell')
        elif int(arcpy.GetCount_management("test_layer").getOutput(0)) == 0:
            logging.error(
                ('\nERROR: ppt_hru_id {0} is not a valid cell hru_id' +
                 '\nERROR: ppt_ratios will not be forced to 1' +
                 ' at cell {0}\n').format(ppt_hru_id))
            ppt_hru_id = 0
        arcpy.Delete_management("test_layer")
        logging.info(
            ('  Observed Mean Monthly PPT ({0}):\n    {1}\n    (Script ' +
             'will assume these are listed in month order, i.e. Jan, ' +
             'Feb, ...)').format(ppt_obs_units, ppt_obs_list))

        # Convert units while reading obs values
        if ppt_obs_units == 'mm':
            factor = 1
        elif ppt_obs_units == 'cm':
            factor = 10
        elif ppt_obs_units == 'm':
            factor = 1000
        elif ppt_obs_units == 'in':
            factor = 25.4
        elif ppt_obs_units == 'ft':
            factor = 304.8
        if factor != 1:
            ppt_obs_list = [p * factor for p in ppt_obs_list]
            logging.info(
                '\n  Converted Mean Monthly PPT ({0}):\n    {1}'.format(
                    ppt_obs_units, ppt_obs_list))

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: Fishnet ({0}) does not exist'.format(
                hru.polygon_path))
        sys.exit()
    if set_ppt_zones_flag:
        if not arcpy.Exists(ppt_zone_orig_path):
            logging.error(
                '\nERROR: PPT Zone ({0}) does not exist'.format(
                    ppt_zone_orig_path))
            sys.exit()
        # ppt_zone_path must be a polygon shapefile
        if arcpy.Describe(ppt_zone_orig_path).datasetType != 'FeatureClass':
            logging.error(
                '\nERROR: ppt_zone_path must be a polygon shapefile')
            sys.exit()
        # Check ppt_zone_fields
        if ppt_zone_field.upper() in ['', 'FID', 'NONE']:
            ppt_zone_field = arcpy.Describe(ppt_zone_orig_path).OIDFieldName
            logging.warning(
                '\n  NOTE: Using {0} to set {1}\n'.format(
                    ppt_zone_field, hru.ppt_zone_id_field))
        elif not arcpy.ListFields(ppt_zone_orig_path, ppt_zone_field):
            logging.error(
                '\nERROR: ppt_zone_field field {0} does not exist\n'.format(
                    ppt_zone_field))
            sys.exit()
        # Need to check that ppt_zone_field is an int type
        elif not [f.type for f in arcpy.Describe(ppt_zone_orig_path).fields
                  if (f.name == ppt_zone_field and
                      f.type in ['SmallInteger', 'Integer'])]:
            logging.error(
                '\nERROR: ppt_zone_field field {0} must be an integer type\n'.format(
                    ppt_zone_field))
            sys.exit()
        # Need to check that ppt_zone_field is all positive values
        elif min([row[0] for row in arcpy.da.SearchCursor(ppt_zone_orig_path, [ppt_zone_field])]) <= 0:
            logging.error(
                '\nERROR: ppt_zone_field values must be positive\n'.format(
                    ppt_zone_field))
            sys.exit()

    # Build output folders if necesssary
    ppt_ratio_temp_ws = os.path.join(hru.param_ws, 'ppt_ratio_temp')
    if not os.path.isdir(ppt_ratio_temp_ws):
        os.mkdir(ppt_ratio_temp_ws)
    ppt_zone_path = os.path.join(ppt_ratio_temp_ws, 'ppt_zone.shp')
    # ppt_zone_clip_path = os.path.join(ppt_ratio_temp_ws, 'ppt_zone_clip.shp')


    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    # env.pyramid = 'PYRAMIDS -1'
    env.pyramid = 'PYRAMIDS 0'
    env.workspace = hru.param_ws
    env.scratchWorkspace = hru.scratch_ws

    # Set month list based on flags
    month_list = ['{0:02d}'.format(m) for m in range(1, 13)]
    ppt_field_list = ['PPT_{0}'.format(m) for m in month_list]
    ratio_field_list = ['PPT_RT_{0}'.format(m) for m in month_list]
    # month_list.extend(['14'])

    # Check fields
    logging.info('\nAdding PRISM fields if necessary')
    # PPT zone fields
    support_functions.add_field_func(
        hru.polygon_path, hru.ppt_zone_id_field, 'LONG')
    # PPT ratio fields
    for ppt_field in ppt_field_list:
        support_functions.add_field_func(hru.polygon_path, ppt_field, 'DOUBLE')

    # Calculate PPT zone ID
    if set_ppt_zones_flag:
        logging.info('\nCalculating cell HRU PPT zone ID')
        ppt_zone_desc = arcpy.Describe(ppt_zone_orig_path)
        ppt_zone_sr = ppt_zone_desc.spatialReference
        logging.debug('  PPT zones: {0}'.format(ppt_zone_orig_path))
        logging.debug('  PPT zones spat. ref.:  {0}'.format(
            ppt_zone_sr.name))
        logging.debug('  PPT zones GCS:         {0}'.format(
            ppt_zone_sr.GCS.name))
        # Reset LAKE_ID
        if set_ppt_zones_flag:
            logging.info('  Resetting {0} to 0'.format(hru.ppt_zone_id_field))
            arcpy.CalculateField_management(
                hru.polygon_path, hru.ppt_zone_id_field, 0, 'PYTHON')
        # If ppt_zone spat_ref doesn't match hru_param spat_ref
        # Project ppt_zone to hru_param spat ref
        # Otherwise, read ppt_zone directly
        if hru.sr.name != ppt_zone_sr.name:
            logging.info('  Projecting PPT zones...')
            # Set preferred transforms
            transform_str = support_functions.transform_func(
                hru.sr, ppt_zone_sr)
            logging.debug('    Transform: {0}'.format(transform_str))
            # Project ppt_zone shapefile
            arcpy.Project_management(
                ppt_zone_orig_path, ppt_zone_path, hru.sr,
                transform_str, ppt_zone_sr)
            del transform_str
        else:
            arcpy.Copy_management(ppt_zone_orig_path, ppt_zone_path)
        # Remove all unnecesary fields
        for field in arcpy.ListFields(ppt_zone_path):
            skip_field_list = ppt_field_list + [ppt_zone_field, 'Shape']
            if field.name not in skip_field_list:
                try:
                    arcpy.DeleteField_management(ppt_zone_path, field.name)
                except:
                    pass
        # Set ppt zone ID
        logging.info('  Setting {0}'.format(hru.ppt_zone_id_field))
        zone_by_centroid_func(
            ppt_zone_path, hru.ppt_zone_id_field, ppt_zone_field,
            hru.polygon_path, hru.point_path, hru_param)
        # zone_by_area_func(
        #    ppt_zone_layer, hru.ppt_zone_id_field, ppt_zone_field,
        #    hru.polygon_path, hru_param, hru_area_field, None, 50)
        # Cleanup
        del ppt_zone_desc, ppt_zone_sr
    else:
        # Set all cells to PPT zone 1
        arcpy.CalculateField_management(
            hru.polygon_path, hru.ppt_zone_id_field, 1, 'PYTHON')


    # Calculate PPT ratios
    logging.info('\nCalculating PRISM mean monthly PPT ratios')
    if set_ppt_zones_flag:
        # Use mean monthly PPT for each zone
        # DEADBEEF - ZONE_VALUE is calculated in zone_by_centroid_func
        # There is probably a cleaner way of linking these two
        fields = ppt_field_list + ['ZONE_VALUE']
        ppt_obs_dict = dict()
        with arcpy.da.SearchCursor(ppt_zone_path, fields) as s_cursor:
            for row in s_cursor:
                # Convert units while reading obs values
                value = map(float, row[:-1])
                if ppt_obs_units == 'mm':
                    value *= 1
                elif ppt_obs_units == 'cm':
                    value *= 10
                elif ppt_obs_units == 'm':
                    value *= 1000
                elif ppt_obs_units == 'in':
                    value *= 25.4
                elif ppt_obs_units == 'ft':
                    value *= 304.8
                ppt_obs_dict[row[-1]] = value
        fields = ppt_field_list + ratio_field_list + [hru.ppt_zone_id_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                for i, month in enumerate(month_list):
                    ppt_i = fields.index('PPT_{0}'.format(month))
                    ratio_i = fields.index('PPT_RT_{0}'.format(month))
                    # If ppt_zone_id was not in zone data, set to 0
                    if row[-1] in ppt_obs_dict.keys():
                        row[ratio_i] = row[ppt_i] / ppt_obs_dict[row[-1]][i]
                    else:
                        row[ratio_i] = 0
                    del ppt_i, ratio_i
                u_cursor.updateRow(row)
            del row
    else:
        # Get PRISM precip at PPT_HRU_ID
        fields = [hru.id_field] + ppt_field_list
        logging.debug('  Fields: {0}'.format(', '.join(fields)))

        # Scale all ratios so PRISM will match observed at a target cell
        if ppt_hru_id != 0:
            ppt_prism_list = map(float, arcpy.da.SearchCursor(
                hru.polygon_path, fields,
                '"{0}" = {1}'.format(hru.id_field, ppt_hru_id)).next()[1:])
            logging.info('  PRISM PPT: {0}'.format(
                ', '.join(['{0:.2f}'.format(p) for p in ppt_prism_list])))
            # Ratio of MEASURED/OBSERVED PPT to PRISM PPT
            # This will be multiplied by PRISM/OBSERVED below
            ppt_ratio_list = [
                float(o) / p if p > 0 else 0
                for o, p in zip(ppt_obs_list, ppt_prism_list)]
            logging.info('  Obs./PRISM: {0}'.format(
                ', '.join(['{0:.3f}'.format(p) for p in ppt_ratio_list])))
        else:
            ppt_ratio_list = [1 for p in ppt_obs_list]

        # Use single mean monthly PPT for all cells
        # Assume ppt_obs_list is in month order
        fields = ppt_field_list + ratio_field_list
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                for i, month in enumerate(month_list):
                    ppt_i = fields.index('PPT_{0}'.format(month))
                    ratio_i = fields.index('PPT_RT_{0}'.format(month))
                    if ppt_obs_list[i] > 0:
                        row[ratio_i] = (
                            ppt_ratio_list[i] * row[ppt_i] / ppt_obs_list[i])
                u_cursor.updateRow(row)
            del row


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='PRISM Precipitation Ratio Parameters',
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

    # Convert relative paths to absolute path
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

    # Calculate GSFLOW PPT Ratio Parameters
    ppt_ratio_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
