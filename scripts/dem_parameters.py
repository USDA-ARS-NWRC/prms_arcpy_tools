#--------------------------------
# Name:         dem_parameters.py
# Purpose:      GSFLOW DEM parameters
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


def dem_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate PRMS DEM Parameters

    Args:
        config_path: Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    # Initialize hru parameters class
    hru = HRUParameters(config_path)

#     # Open input parameter config file
#     inputs_cfg = ConfigParser.ConfigParser()
#     try:
#         inputs_cfg.readfp(open(config_path))
#     except:
#         logging.error('\nERROR: Config file could not be read, ' +
#                       'is not an input file, or does not exist\n' +
#                       'ERROR: config_file = {0}\n').format(config_path)
#         sys.exit()
#     logging.debug('\nReading Input File')

    # Log DEBUG to file
    log_file_name = 'dem_parameters_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nPRMS DEM Parameters')

    # read the DEM parameters from the file, check if DEM exists
    hru.read_DEM_parameters()
    
    if hru.calc_flow_acc_dem_flag:
        # Get factor for scaling dem_flowacc values to avoid 32 bit int limits
        try:
            flow_acc_dem_factor = float(hru.inputs_cfg.get('INPUTS', 'flow_acc_dem_factor'))
        except:
            # This is a worst case for keeping flow_acc_dem from exceeding 2E9
            # Assume all cells flow to 1 cell
            flow_acc_dem_factor = int(
                arcpy.GetCount_management(hru.point_path).getOutput(0))
            # Assume flow acc is in every DEM cell in HRU cell
            flow_acc_dem_factor *= (float(300.0) / hru.dem_cs) ** 2
            # Need to account for the elevation in this worst cell
            # For now just make it 100
            # flow_acc_dem_factor *= max_elevation
            flow_acc_dem_factor *= 100
            # Calculate ratio of flow_acc_dem to a 32 bit int
            flow_acc_dem_factor /= (0.5 * 2**32)
            # If the ratio is less than 0.1, round up to 0.1 so factor -> 1.0
            flow_acc_dem_factor = min(0.1, flow_acc_dem_factor)
            # Round up to next multiple of 10 just to be safe
            flow_acc_dem_factor = 1.0 / 10**(int(math.log10(flow_acc_dem_factor))+1)
            logging.info(
                ('  flow_acc_dem_factor was not set in the input file\n' +
                 '  Using automatic flow_acc_dem_factor: {0}').format(
                     flow_acc_dem_factor))

    # Calc flow_acc/flow_dir
    # DEADBEEF - For now, set these to True only if needed
    # calc_flow_acc_flag = inputs_cfg.getboolean('INPUTS', 'calc_flow_acc_flag')
    # calc_flow_dir_flag = inputs_cfg.getboolean('INPUTS', 'calc_flow_dir_flag')
    if hru.calc_flow_acc_dem_flag:
        calc_flow_acc_flag = True
        calc_flow_dir_flag = True
    else:
        calc_flow_acc_flag = False
        calc_flow_dir_flag = False

    # Check input paths
    hru.check_polygon_path()

    # Read remap from config file and check that remap folder is valid
    hru.read_remap_parameters()

    # Check that remap files exist
    # Check remap files comment style
    aspect_remap_path = os.path.join(hru.remap_ws, hru.aspect_remap_name)
    temp_adj_remap_path = os.path.join(hru.remap_ws, hru.temp_adj_remap_name)
    remap_path_list = [aspect_remap_path, temp_adj_remap_path]
    for remap_path in remap_path_list:
        remap_check(remap_path)

    # DEADBEEF
    # if not os.path.isfile(aspect_remap_path):
    #    logging.error(
    #        '\nERROR: ASCII remap file ({0}) does not exist\n'.format(
    #            os.path.basename(aspect_remap_path)))
    #    sys.exit()
    # if not os.path.isfile(temp_adj_remap_path):
    #    logging.error(
    #        '\nERROR: ASCII remap file ({0}) does not exist\n'.format(
    #            os.path.basename(temp_adj_remap_path)))
    #    sys.exit()
    #  Check remap files comment style
    # if '10.2' in arcpy.GetInstallInfo()['version']:
    #    if remap_comment_check(aspect_remap_path):
    #        logging.error(
    #            ('\nERROR: ASCII remap file ({0}) has pre-ArcGIS 10.2 ' +
    #             'comments\n').format(os.path.basename(aspect_remap_path)))
    #        sys.exit()
    #    if remap_comment_check(temp_adj_remap_path):
    #        logging.error(
    #            ('\nERROR: ASCII remap file ({0}) has pre-ArcGIS 10.2 ' +
    #             'comments\n').format(os.path.basename(temp_adj_remap_path)))
    #        sys.exit()


    # Build output folder if necessary
    dem_temp_ws = os.path.join(hru.param_ws, 'dem_rasters')
    if not os.path.isdir(dem_temp_ws):
        os.mkdir(dem_temp_ws)

    # Output paths
    dem_path = os.path.join(dem_temp_ws, 'dem.img')
    dem_fill_path = os.path.join(dem_temp_ws, 'dem_fill.img')
    flow_dir_path = os.path.join(dem_temp_ws, 'flow_dir.img')
    flow_acc_path = os.path.join(dem_temp_ws, 'flow_acc.img')
    flow_acc_dem_path = os.path.join(dem_temp_ws, 'flow_acc_x_dem.img')
    flow_acc_filter_path = os.path.join(dem_temp_ws, 'flow_acc_filter.img')
    dem_integer_path = os.path.join(dem_temp_ws, 'dem_integer.img')
    dem_slope_path = os.path.join(dem_temp_ws, 'dem_slope.img')
    dem_aspect_path = os.path.join(dem_temp_ws, 'dem_aspect.img')
    dem_aspect_reclass_path = os.path.join(dem_temp_ws, 'aspect_reclass.img')
    temp_adj_path = os.path.join(dem_temp_ws, 'temp_adj.img')

    # Set ArcGIS environment variables
    arcpy.CheckOutExtension('Spatial')
    env.overwriteOutput = True
    env.pyramid = 'PYRAMIDS -1'
    # env.pyramid = 'PYRAMIDS 0'
    # env.rasterStatistics = 'NONE'
    # env.extent = 'MINOF'
    env.workspace = dem_temp_ws
    env.scratchWorkspace = hru.scratch_ws

    # Check DEM field
    logging.info('\nAdding DEM fields if necessary')
    add_field_func(hru.polygon_path, hru.dem_mean_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.dem_median_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.dem_max_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.dem_min_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.dem_adj_field, 'DOUBLE')
    if hru.calc_flow_acc_dem_flag:
        add_field_func(hru.polygon_path, hru.dem_flowacc_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_sum_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_count_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.dem_sink8_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.dem_sink4_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.elev_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.aspect_field, 'LONG')
    add_field_func(hru.polygon_path, hru.slope_deg_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.slope_rad_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.slope_pct_field, 'DOUBLE')
    # add_field_func(hru.polygon_path, hru.deplcrv_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.jh_tmin_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.jh_tmax_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.jh_coef_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.snarea_thresh_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.tmax_adj_field, 'DOUBLE')
    add_field_func(hru.polygon_path, hru.tmin_adj_field, 'DOUBLE')

    

    # Assume all DEM rasters will need to be rebuilt
    # Check slope, aspect, and proejcted DEM rasters
    # This will check for matching spat. ref., snap point, and cellsize

    # If DEM is GCS, project it to 10m to match
    # DEADBEEF - I had originally wanted the DEM to get projected only once
    #   but if the user wants to rerun this script, then all steps should
    #   be rerun.  This also allows the user to change the DEM raster
    # dem_flag = valid_raster_func(
    #    dem_path, 'projected DEM', hru, dem_cs)
    # if arcpy.Exists(dem_orig_path) and not dem_flag:
    logging.info('\nProjecting DEM raster')
    dem_orig_sr = Raster(hru.dem_orig_path).spatialReference
    logging.debug('  DEM GCS:   {0}'.format(
        dem_orig_sr.GCS.name))
    
    # Remove existing projected DEM
    if arcpy.Exists(dem_path):
        arcpy.Delete_management(dem_path)
        
    # Set preferred transforms
    transform_str = transform_func(hru.sr, dem_orig_sr)
    logging.debug('  Transform: {0}'.format(transform_str))
    logging.debug('  Projection method: {0}'.format(hru.dem_proj_method))
    
    # Project DEM
    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    logging.debug('  Input SR:  {0}'.format(dem_orig_sr.exportToString()))
    logging.debug('  Output SR: {0}'.format(hru.sr.exportToString()))
    project_raster_func(
        hru.dem_orig_path, dem_path, hru.sr,
        hru.dem_proj_method, hru.dem_cs, transform_str,
        dem_orig_sr, hru)
    # env.extent = hru.extent
    # arcpy.ProjectRaster_management(
    #    dem_orig_path, dem_path, hru.sr,
    #    dem_proj_method, dem_cs, transform_str,
    #    '{0} {1}'.format(hru.ref_x, hru.ref_y),
    #    dem_orig_sr)
    # arcpy.ClearEnvironment('extent')

    # Check linear unit of raster
    # DEADBEEF - The conversion could probably be dynamic
    dem_obj = Raster(dem_path)
    linear_unit_list = ['METER', 'FOOT_US', 'FOOT']
    linear_unit = dem_obj.spatialReference.linearUnitName.upper()
    if linear_unit not in linear_unit_list:
        logging.error(
            '\nERROR: The linear unit of the projected/clipped DEM must' +
            ' be meters or feet\n  {0}'.format(linear_unit))
        sys.exit()
    del dem_obj

    # Calculate filled DEM, flow_dir, & flow_acc
    logging.info('\nCalculating filled DEM raster')
    dem_fill_obj = Fill(dem_path)
    dem_fill_obj.save(dem_fill_path)
    del dem_fill_obj
    if calc_flow_dir_flag:
        logging.info('Calculating flow direction raster')
        dem_fill_obj = Raster(dem_fill_path)
        flow_dir_obj = FlowDirection(dem_fill_obj, True)
        flow_dir_obj.save(flow_dir_path)
        del flow_dir_obj, dem_fill_obj
    if calc_flow_acc_flag:
        logging.info('Calculating flow accumulation raster')
        flow_dir_obj = Raster(flow_dir_path)
        flow_acc_obj = FlowAccumulation(flow_dir_obj)
        flow_acc_obj.save(flow_acc_path)
        del flow_acc_obj, flow_dir_obj
    if hru.calc_flow_acc_dem_flag:
        # flow_acc_dem_obj = dem_fill_obj * flow_acc_obj
        # Low pass filter of flow_acc then take log10
        flow_acc_filter_obj = Filter(Raster(flow_acc_path), 'LOW', 'NODATA')
        flow_acc_filter_obj *= flow_acc_dem_factor
        flow_acc_filter_obj.save(flow_acc_filter_path)
        flow_acc_dem_obj = Raster(dem_fill_path) * flow_acc_filter_obj
        flow_acc_dem_obj.save(flow_acc_dem_path)
        del flow_acc_dem_obj, flow_acc_filter_obj

    # Calculate an integer version of DEM for median zonal stats
    dem_integer_obj = Int(Raster(dem_path) * 100)
    dem_integer_obj.save(dem_integer_path)
    del dem_integer_obj

    # Calculate slope
    logging.info('Calculating slope raster')
    dem_slope_obj = Slope(dem_fill_path, 'DEGREE')
    # Setting small slopes to zero
    logging.info('  Setting slopes <= 0.01 to 0')
    dem_slope_obj = Con(dem_slope_obj <= 0.01, 0, dem_slope_obj)
    dem_slope_obj.save(dem_slope_path)
    del dem_slope_obj

    # Calculate aspect
    logging.info('Calculating aspect raster')
    dem_aspect_obj = Aspect(dem_fill_path)
    # Set small slopes to -1 aspect
    logging.debug('  Setting aspect for slopes <= 0.01 to -1')
    dem_aspect_obj = Con(Raster(dem_slope_path) > 0.01, dem_aspect_obj, -1)
    dem_aspect_obj.save(dem_aspect_path)
    del dem_aspect_obj

    # Reclassify aspect
    logging.debug('  Reclassifying: {0}'.format(aspect_remap_path))
    dem_aspect_reclass_obj = ReclassByASCIIFile(
        dem_aspect_path, aspect_remap_path)
    dem_aspect_reclass_obj.save(dem_aspect_reclass_path)
    del dem_aspect_reclass_obj

    # Temperature Aspect Adjustment
    logging.info('Calculating temperature aspect adjustment raster')
    temp_adj_obj = Float(ReclassByASCIIFile(
        dem_aspect_reclass_path, temp_adj_remap_path))
    # Since reclass can't remap to floats directly
    # Values are scaled by 10 and stored as integers
    temp_adj_obj *= 0.1
    # This is a function in the
    # temp_adj_obj = reclass_ascii_float_func(
    #    dem_aspect_reclass_path, temp_adj_remap_path)
    temp_adj_obj.save(temp_adj_path)
    del temp_adj_obj


    # List of rasters, fields, and stats for zonal statistics
    zs_dem_dict = dict()
    zs_dem_dict[hru.dem_mean_field] = [dem_path, 'MEAN']
    if hru.calc_flow_acc_dem_flag:
        zs_dem_dict[hru.dem_sum_field] = [flow_acc_dem_path, 'SUM']
        zs_dem_dict[hru.dem_count_field] = [flow_acc_filter_path, 'SUM']
    # CGM - Zonal stats wasn't working with median
    # zs_dem_dict[hru.dem_median_field] = [dem_integer_path, 'MEDIAN']
    zs_dem_dict[hru.dem_max_field] = [dem_path, 'MAXIMUM']
    zs_dem_dict[hru.dem_min_field] = [dem_path, 'MINIMUM']
    # zs_dem_dict[hru.elev_field]   = [dem_integer_path, 'MEDIAN']
    # zs_dem_dict[hru.aspect_field] = [dem_aspect_path, 'MINIMUM']
    zs_dem_dict[hru.aspect_field] = [dem_aspect_reclass_path, 'MAJORITY']
    zs_dem_dict[hru.slope_deg_field] = [dem_slope_path, 'MEAN']
    zs_dem_dict[hru.tmax_adj_field] = [temp_adj_path, 'MEAN']
    zs_dem_dict[hru.tmin_adj_field] = [temp_adj_path, 'MEAN']


    # Calculate DEM zonal statistics
    logging.info('\nCalculating DEM zonal statistics')
    zonal_stats_func(zs_dem_dict, hru.polygon_path, hru.point_path, hru)

    # Reset DEM_MEDIAN
    # logging.info('\nCalculating {0}'.format(hru.dem_median_field))
    # arcpy.CalculateField_management(
    #    hru.polygon_path, hru.dem_median_field,
    #    # Convert meters to feet
    #    '0.01 * !{0}!'.format(hru.dem_median_field), 'PYTHON')

    # Calculate HRU_ELEV (HRU elevation in feet)
    logging.info('\nCalculating initial {0} from {1}'.format(
        hru.elev_field, hru.dem_mean_field))
    if linear_unit in ['METERS']:
        logging.info('  Converting from meters to feet')
        arcpy.CalculateField_management(
            hru.polygon_path, hru.elev_field,
            '!{0}! * 3.28084'.format(hru.dem_mean_field), 'PYTHON')
    elif linear_unit in ['FOOT_US', 'FOOT']:
        arcpy.CalculateField_management(
            hru.polygon_path, hru.elev_field,
            '!{0}!'.format(hru.dem_mean_field), 'PYTHON')


    # Flow accumulation weighted elevation
    if hru.calc_flow_acc_dem_flag:
        logging.info('Calculating {0}'.format(hru.dem_flowacc_field))
        hru_polygon_layer = 'hru_polygon_layer'
        arcpy.MakeFeatureLayer_management(
            hru.polygon_path, hru_polygon_layer)
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" > 0'.format(hru.dem_count_field))
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.dem_flowacc_field,
            'float(!{0}!) / !{1}!'.format(hru.dem_sum_field, hru.dem_count_field),
            'PYTHON')
        # Clear dem_flowacc for any cells that have zero sum or count
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '("{0}" = 0) OR ("{1}" = 0)'.format(
                hru.dem_count_field, hru.dem_sum_field))
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.dem_flowacc_field, 0, 'PYTHON')
        arcpy.Delete_management(hru_polygon_layer)
        # arcpy.DeleteField_management(hru.polygon_path, hru.dem_sum_field)
        # arcpy.DeleteField_management(hru.polygon_path, hru.dem_count_field)

    # Fill DEM_ADJ if it is not set
    if all([row[0] == 0 for row in arcpy.da.SearchCursor(
        hru.polygon_path, [hru.dem_adj_field])]):
        logging.info('Filling {0} from {1}'.format(
            hru.dem_adj_field, hru.dem_adj_copy_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.dem_adj_field,
            'float(!{0}!)'.format(hru.dem_adj_copy_field), 'PYTHON')
    elif hru.reset_dem_adj_flag:
        logging.info('Filling {0} from {1}'.format(
            hru.dem_adj_field, hru.dem_adj_copy_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.dem_adj_field,
            'float(!{0}!)'.format(hru.dem_adj_copy_field), 'PYTHON')
    else:
        logging.info(
            ('{0} appears to already have been set and ' +
             'will not be overwritten').format(hru.dem_adj_field))

    # HRU_SLOPE in radians
    logging.info('Calculating {0} (Slope in Radians)'.format(
        hru.slope_rad_field))
    arcpy.CalculateField_management(
        hru.polygon_path, hru.slope_rad_field,
        'math.pi * !{0}! / 180'.format(hru.slope_deg_field), 'PYTHON')
    # HRU_SLOPE in percent
    logging.info('Calculating {0} (Percent Slope)'.format(
        hru.slope_pct_field))
    arcpy.CalculateField_management(
        hru.polygon_path, hru.slope_pct_field,
        'math.tan(!{0}!)'.format(hru.slope_rad_field), 'PYTHON')

    # HRU_DEPLCRV
    # deplcrv is set to 1 for all active cells when writing parameter file
    # logging.info('Calculating {0}'.format(hru.deplcrv_field))
    # arcpy.CalculateField_management(
    #    hru.polygon_path, hru.deplcrv_field, '1', 'PYTHON')

    # Jensen-Haise Potential ET air temperature coefficient
    logging.info('Calculating JH_COEF_HRU')
    # First check if PRISM TMAX/TMIN have been set
    # If max July value is 0, use default values
    if (hru.calc_prism_jh_coef_flag and
        (len(arcpy.ListFields(hru.polygon_path, 'TMAX_07')) == 0 or
         field_stat_func(hru.polygon_path, 'TMAX_07', 'MAXIMUM') == 0)):
        hru.calc_prism_jh_coef_flag = False
        
    # Use PRISM temperature values
    if hru.calc_prism_jh_coef_flag:
        logging.info('  Using PRISM temperature values')
        tmax_field_list = ['!TMAX_{0:02d}!'.format(m) for m in range(1, 13)]
        tmin_field_list = ['!TMIN_{0:02d}!'.format(m) for m in range(1, 13)]
        tmax_expr = 'max([{0}])'.format(','.join(tmax_field_list))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.jh_tmax_field, tmax_expr, 'PYTHON')
        # Get TMIN for same month as maximum TMAX
        tmin_expr = 'max(zip([{0}],[{1}]))[1]'.format(
            ','.join(tmax_field_list), ','.join(tmin_field_list))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.jh_tmin_field, tmin_expr, 'PYTHON')
    # Use default temperature values
    else:
        logging.info('  Using default temperature values (7 & 25)')
        arcpy.CalculateField_management(
            hru.polygon_path, hru.jh_tmax_field, 25, 'PYTHON')
        arcpy.CalculateField_management(
            hru.polygon_path, hru.jh_tmin_field, 7, 'PYTHON')
    jensen_haise_func(
        hru.polygon_path, hru.jh_coef_field, hru.elev_field,
        hru.jh_tmin_field, hru.jh_tmax_field)

    # SNAREA_THRESH
    logging.info('Calculating {0}'.format(hru.snarea_thresh_field))
    elev_min = field_stat_func(hru.polygon_path, hru.elev_field, 'MINIMUM')
    arcpy.CalculateField_management(
        hru.polygon_path, hru.snarea_thresh_field,
        '(!{0}! - {1}) * 0.005'.format(hru.elev_field, elev_min),
        'PYTHON')

    # Clear slope/aspect values for lake cells (HRU_TYPE == 2)
    # Also clear for ocean cells (HRU_TYPE == 0 and DEM_ADJ == 0)
    if True:
        logging.info('\nClearing slope/aspect parameters for lake cells')
        hru_polygon_layer = "hru_polygon_layer"
        arcpy.MakeFeatureLayer_management(
            hru.polygon_path, hru_polygon_layer)
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 2 OR ("{0}" = 0 AND "{1}" = 0)'.format(
                hru.type_field, hru.dem_adj_field))
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.aspect_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.slope_deg_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.slope_rad_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.slope_pct_field, 0, 'PYTHON')
        # arcpy.CalculateField_management(
        #    hru_polygon_layer, hru.deplcrv_field, 0, 'PYTHON')
        # arcpy.CalculateField_management(
        #    hru_polygon_layer, hru.snarea_field, 0, 'PYTHON')
        # arcpy.CalculateField_management(
        #    hru_polygon_layer, hru.tmax_adj_field, 0, 'PYTHON')
        # arcpy.CalculateField_management(
        #    hru_polygon_layer, hru.tmin_adj_field, 0, 'PYTHON')

        # Should JH coefficients be cleared for lakes?
        # logging.info('\nClearing JH parameters for ocean cells')
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 0 AND "{1}" = 0'.format(
                hru.type_field, hru.dem_adj_field))
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.jh_coef_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.jh_tmax_field, 0, 'PYTHON')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.jh_tmin_field, 0, 'PYTHON')

        arcpy.Delete_management(hru_polygon_layer)
        del hru_polygon_layer
        
    logging.info('Done!')


def field_stat_func(input_path, value_field, stat='MAXIMUM'):
    """"""
    value_list = []
    with arcpy.da.SearchCursor(input_path, value_field) as s_cursor:
        for row in s_cursor:
            value_list.append(row[0])
    if stat.upper() in ['MAXIMUM', 'MAX']:
        return max(value_list)
    elif stat.upper() in ['MINIMUM', 'MIN']:
        return min(value_list)
    else:
        return float(sum(value_list)) / count(value_list)


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='DEM Parameters',
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

    # Calculate GSFLOW DEM Parameters
    dem_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
