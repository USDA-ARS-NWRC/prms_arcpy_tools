#--------------------------------
# Name:         support_functions.py
# Purpose:      GSFLOW parameter support functions
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2016-02-26
# Python:       2.7
#--------------------------------

from collections import defaultdict
import ConfigParser
import heapq
import itertools
import logging
import math
from operator import itemgetter
import os
import re
import sys
from time import sleep

import numpy as np

import arcpy
from arcpy import env
from arcpy.sa import *

from support_functions import *


class HRUParameters():
    """"""
    def __init__(self, config_path):
        
        # Open input parameter config file
        inputs_cfg = ConfigParser.ConfigParser()
        try:
            inputs_cfg.readfp(open(config_path))
        except IOError:
            logging.error(('\nERROR: Config file does not exist\n'))
            sys.exit()
        except ConfigParser.MissingSectionHeaderError:
            logging.error('\nERROR: Config file is missing a section header\n' +
                          '    Please make sure the following line is at the ' +
                          'beginning of the file\n[INPUTS]\n')
            sys.exit()
        except:
            logging.error(('\nERROR: Config file could not be read\n' +
                           '  {0}\n').format(config_path))
        logging.debug('\nReading Input File')
        logging.debug('  {}'.format(os.path.basename(config_path)))

        # Open field list config file
        # Use script directory (from sys.argv[0]) in case script is a
        #   relative path (i.e. called from a project folder)
        field_list_path = os.path.join(
            os.path.dirname(sys.argv[0]), 'field_list.ini')
        # #field_list_path =  inputs_cfg.get('INPUTS', 'field_list_path')
        fields_cfg = ConfigParser.ConfigParser()
        try:
            fields_cfg.readfp(open(field_list_path))
        except IOError:
            logging.error(('\nERROR: Field list file does not exist\n' +
                           '  {0}\n').format(field_list_path))
            sys.exit()
        except ConfigParser.MissingSectionHeaderError:
            logging.error('\nERROR: Field list file is missing a section header\n' +
                          '    Please make sure the following line is at the ' +
                          'beginning of the file\n[FIELDS]\n')
            sys.exit()
        except:
            logging.error(('\nERROR: Field list file could not be read\n' +
                           '  {0}\n').format(field_list_path))
        logging.debug('\nReading Field List File')

        # Read parameters from config file
        self.polygon_path = inputs_cfg.get('INPUTS', 'hru_path')
        self.point_path   = inputs_cfg.get('INPUTS', 'hru_centroid_path')
        self.sr_name = inputs_cfg.get('INPUTS', 'hru_projection')
        self.snap_method = inputs_cfg.get('INPUTS', 'hru_param_snap_method')
        self.fid_field = inputs_cfg.get('INPUTS', 'orig_fid_field')
        self.type_field = fields_cfg.get('FIELDS', 'type_field')

        self.param_ws = inputs_cfg.get('INPUTS', 'parameter_folder')
        if not os.path.isdir(self.param_ws):
            os.mkdir(self.param_ws)

        # Log workspace
        self.log_ws = os.path.join(self.param_ws, 'logs')
        if not os.path.isdir(self.log_ws):
            os.mkdir(self.log_ws)

        # Scratch workspace
        try:
            scratch_name = inputs_cfg.get('INPUTS', 'scratch_name')
        except:
            scratch_name = 'in_memory'
        if scratch_name == 'in_memory':
            self.scratch_ws = scratch_name
        else:
            scratch_ws = os.path.join(self.param_ws, scratch_name)
            if not os.path.isdir(scratch_ws):
                os.mkdir(scratch_ws)
            self.scratch_ws = scratch_ws

        # Set spatial reference of hru shapefile
        if arcpy.Exists(self.polygon_path):
            hru_desc = arcpy.Describe(self.polygon_path)
            self.sr = hru_desc.spatialReference
            self.extent = round_extent(hru_desc.extent, 6)
            logging.info('  HRU extent:     {0}'.format(extent_string(self.extent)))
            logging.debug('  HRU spat. ref.: {0}'.format(self.sr.name))
            logging.debug('  HRU GCS:        {0}'.format(self.sr.GCS.name))

#             # Check that the HRU is snapped to the reference point
#             if not snapped(self.extent, self.ref_pnt, self.cs):
#                 logging.error(
#                     ('\nWARNING: {0} does not appear to be snapped to the INI ' +
#                      'file reference point\n  This may be a rounding issue.').format(
#                         os.path.basename(self.polygon_path)))
# #                 raw_input('Press ENTER to continue')

            # DEADBEEF - I'm not sure why I would adjust the extent
            # If the extent doesn't match the refence point, the script
            #   should probably terminate
            # self.extent = adjust_extent_to_snap(
            #    hru_desc.extent, self.ref_pnt, self.cs, 'ROUND')
            # self.extent = adjust_extent_to_snap(
            #    hru_param_desc.extent, snap_pnt, self.cs, 'ROUND')

        # Some fields are dependent on the control flags
        try:
            self.set_lake_flag = inputs_cfg.getboolean('INPUTS', 'set_lake_flag')
        except:
            logging.debug('  set_lake_flag = False')
            self.set_lake_flag = False
        
        self.calc_flow_acc_dem_flag = inputs_cfg.getboolean('INPUTS', 'calc_flow_acc_dem_flag')
        self.calc_topo_index_flag = inputs_cfg.getboolean('INPUTS', 'calc_topo_index_flag')
        self.clip_root_depth_flag = inputs_cfg.getboolean('INPUTS', 'clip_root_depth_flag')
        # set_ppt_zones_flag = inputs_cfg.getboolean('INPUTS', 'set_ppt_zones_flag')
#         self.calc_layer_thickness_flag = inputs_cfg.getboolean(
#             'INPUTS', 'calc_layer_thickness_flag')

        # Read in all field names
        self.id_field = fields_cfg.get('FIELDS', 'id_field')
        self.type_in_field = fields_cfg.get('FIELDS', 'type_in_field')
        self.type_field = fields_cfg.get('FIELDS', 'type_field')
        self.dem_mean_field = fields_cfg.get('FIELDS', 'dem_mean_field')
        self.dem_median_field = fields_cfg.get('FIELDS', 'dem_median_field')
        self.dem_max_field = fields_cfg.get('FIELDS', 'dem_max_field')
        self.dem_min_field = fields_cfg.get('FIELDS', 'dem_min_field')
        self.dem_adj_field = fields_cfg.get('FIELDS', 'dem_adj_field')
        
        if self.calc_flow_acc_dem_flag:
            # self.dem_sum_field = 'DEM_SUM'
            # self.dem_count_field = 'DEM_COUNT'
            self.dem_sum_field = fields_cfg.get('FIELDS', 'dem_sum_field')
            self.dem_count_field = fields_cfg.get('FIELDS', 'dem_count_field')
            self.dem_flowacc_field = fields_cfg.get('FIELDS', 'dem_flowacc_field')
        else:
            self.dem_sum_field = 'DEM_SUM'
            self.dem_count_field = 'DEM_COUNT'
            self.dem_flowacc_field = 'DEM_FLOW_AC'
        
        self.dem_sink8_field = fields_cfg.get('FIELDS', 'dem_sink8_field')
        self.dem_sink4_field = fields_cfg.get('FIELDS', 'dem_sink4_field')
        self.area_field = fields_cfg.get('FIELDS', 'area_field')
        self.elev_field = fields_cfg.get('FIELDS', 'elev_field')
        self.aspect_field = fields_cfg.get('FIELDS', 'aspect_field')
        self.slope_deg_field = fields_cfg.get('FIELDS', 'slope_deg_field')
        self.slope_rad_field = fields_cfg.get('FIELDS', 'slope_rad_field')
        self.slope_pct_field = fields_cfg.get('FIELDS', 'slope_pct_field')
        self.topo_index_field = fields_cfg.get('FIELDS', 'topo_index_field')
        self.x_field = fields_cfg.get('FIELDS', 'x_field')
        self.y_field = fields_cfg.get('FIELDS', 'y_field')
        self.lat_field = fields_cfg.get('FIELDS', 'lat_field')
        self.lon_field = fields_cfg.get('FIELDS', 'lon_field')
        self.xlong_field = fields_cfg.get('FIELDS', 'xlong_field')
        self.ylat_field = fields_cfg.get('FIELDS', 'ylat_field')

        if self.set_lake_flag:
            self.lake_id_field = fields_cfg.get('FIELDS', 'lake_id_field')
            self.lake_area_field = fields_cfg.get('FIELDS', 'lake_area_field')
        else:
            self.lake_id_field = 'LAKE_ID'
            self.lake_area_field = 'LAKE_AREA'

        # DEM based
        # self.deplcrv_field = fields_cfg.get('FIELDS', 'deplcrv_field')
        self.jh_tmax_field = fields_cfg.get('FIELDS', 'jh_tmax_field')
        self.jh_tmin_field = fields_cfg.get('FIELDS', 'jh_tmin_field')
        self.jh_coef_field = fields_cfg.get('FIELDS', 'jh_coef_field')
        self.snarea_thresh_field = fields_cfg.get('FIELDS', 'snarea_thresh_field')
        self.tmax_adj_field = fields_cfg.get('FIELDS', 'tmax_adj_field')
        self.tmin_adj_field = fields_cfg.get('FIELDS', 'tmin_adj_field')

        # Vegetation
        self.cov_type_field = fields_cfg.get('FIELDS', 'cov_type_field')
        self.covden_sum_field = fields_cfg.get('FIELDS', 'covden_sum_field')
        self.covden_win_field = fields_cfg.get('FIELDS', 'covden_win_field')
        self.snow_intcp_field = fields_cfg.get('FIELDS', 'snow_intcp_field')
        self.wrain_intcp_field = fields_cfg.get('FIELDS', 'wrain_intcp_field')
        self.srain_intcp_field = fields_cfg.get('FIELDS', 'srain_intcp_field')
        self.rad_trncf_field = fields_cfg.get('FIELDS', 'rad_trncf_field')

        # Soil
        self.awc_field = fields_cfg.get('FIELDS', 'awc_field')
        self.clay_pct_field = fields_cfg.get('FIELDS', 'clay_pct_field')
        self.sand_pct_field = fields_cfg.get('FIELDS', 'sand_pct_field')
        # self.silt_pct_field = fields_cfg.get('FIELDS', 'silt_pct_field')
        self.ksat_field = fields_cfg.get('FIELDS', 'ksat_field')
        self.soil_depth_field = fields_cfg.get('FIELDS', 'soil_depth_field')
        self.root_depth_field = fields_cfg.get('FIELDS', 'root_depth_field')
        self.soil_type_field  = fields_cfg.get('FIELDS', 'soil_type_field')
        self.moist_init_field = fields_cfg.get('FIELDS', 'moist_init_field')
        self.moist_max_field  = fields_cfg.get('FIELDS', 'moist_max_field')
        self.rechr_init_field = fields_cfg.get('FIELDS', 'rechr_init_field')
        self.rechr_max_field  = fields_cfg.get('FIELDS', 'rechr_max_field')
        self.ssr2gw_rate_field = fields_cfg.get('FIELDS', 'ssr2gw_rate_field')
        # self.pref_flow_den_field = fields_cfg.get('FIELDS', 'pref_flow_den_field')
        self.slowcoef_lin_field = fields_cfg.get('FIELDS', 'slowcoef_lin_field')
        self.slowcoef_sq_field  = fields_cfg.get('FIELDS', 'slowcoef_sq_field')
        self.fastcoef_lin_field = fields_cfg.get('FIELDS', 'fastcoef_lin_field')
        self.fastcoef_sq_field  = fields_cfg.get('FIELDS', 'fastcoef_sq_field')

        # Impervious Parameter Fields
        self.imperv_pct_field = fields_cfg.get('FIELDS', 'imperv_pct_field')
        # carea_min_field  = fields_cfg.get('FIELDS', 'carea_min_field')
        self.carea_max_field  = fields_cfg.get('FIELDS', 'carea_max_field')

        # Streams
        self.irunbound_field  = fields_cfg.get('FIELDS', 'irunbound_field')
        self.iseg_field       = fields_cfg.get('FIELDS', 'iseg_field')
        self.flow_dir_field   = fields_cfg.get('FIELDS', 'flow_dir_field')
        self.krch_field       = fields_cfg.get('FIELDS', 'krch_field')
        self.irch_field       = fields_cfg.get('FIELDS', 'irch_field')
        self.jrch_field       = fields_cfg.get('FIELDS', 'jrch_field')
        self.reach_field      = fields_cfg.get('FIELDS', 'reach_field')
        self.rchlen_field     = fields_cfg.get('FIELDS', 'rchlen_field')
        self.maxreach_field   = fields_cfg.get('FIELDS', 'maxreach_field')
        self.outseg_field     = fields_cfg.get('FIELDS', 'outseg_field')
        self.iupseg_field     = fields_cfg.get('FIELDS', 'iupseg_field')
        self.strm_top_field   = fields_cfg.get('FIELDS', 'strm_top_field')
        self.strm_slope_field = fields_cfg.get('FIELDS', 'strm_slope_field')
        self.subbasin_field   = fields_cfg.get('FIELDS', 'subbasin_field')
        self.segbasin_field   = fields_cfg.get('FIELDS', 'segbasin_field')
        self.outflow_field    = fields_cfg.get('FIELDS', 'outflow_field')
        self.hru_segment      = fields_cfg.get('FIELDS', 'hru_segment')
        self.k_coef           = fields_cfg.get('FIELDS', 'k_coef')
        self.obsin_segment    = fields_cfg.get('FIELDS', 'obsin_segment')
        self.tosegment        = fields_cfg.get('FIELDS', 'tosegment')
        self.x_coef           = fields_cfg.get('FIELDS', 'x_coef')
        self.stream_path      = inputs_cfg.get('INPUTS', 'streams_path')


        # if set_ppt_zones_flag:
        self.ppt_zone_id_field = fields_cfg.get('FIELDS', 'ppt_zone_id_field')

#         # Calculate layer thickness and bottoms
#         if self.calc_layer_thickness_flag:
#             self.alluv_field = fields_cfg.get('FIELDS', 'alluv_field')
#             self.alluv_thick_field = fields_cfg.get('FIELDS', 'alluv_thick_field')
#             self.lay1_thick_field = fields_cfg.get('FIELDS', 'lay1_thick_field')
#             self.lay2_thick_field = fields_cfg.get('FIELDS', 'lay2_thick_field')
#             self.lay3_thick_field = fields_cfg.get('FIELDS', 'lay3_thick_field')
#             self.lay4_thick_field = fields_cfg.get('FIELDS', 'lay4_thick_field')
#             self.lay1_bottom_field = fields_cfg.get('FIELDS', 'lay1_bottom_field')
#             self.lay2_bottom_field = fields_cfg.get('FIELDS', 'lay2_bottom_field')
#             self.lay3_bottom_field = fields_cfg.get('FIELDS', 'lay3_bottom_field')
#             self.lay4_bottom_field = fields_cfg.get('FIELDS', 'lay4_bottom_field')

        
        # save the config file
        self.inputs_cfg = inputs_cfg
        
    def check_polygon_path(self):
        """
        Check that the polygon path exists 
        """
        if not arcpy.Exists(self.polygon_path):
            logging.error(
                '\nERROR: HRU ({0}) does not exist\n'.format(self.polygon_path))
            sys.exit()
            
    def read_remap_parameters(self):
        """
        Read the remap parameters from the config file and ensure
        that the folder exists
        """
        
        # Remap
        self.remap_ws = self.inputs_cfg.get('INPUTS', 'remap_folder')
        try:
            self.aspect_remap_name = self.inputs_cfg.get('INPUTS', 'aspect_remap')
        except:
            self.aspect_remap_name = None
            
        try:
            self.temp_adj_remap_name = self.inputs_cfg.get('INPUTS', 'temp_adj_remap')
        except:
            self.temp_adj_remap_name = None
            
        try:
            self.cov_type_remap_name = self.inputs_cfg.get('INPUTS', 'cov_type_remap')
        except:
            self.cov_type_remap_name = None
            
        try:
            self.covden_sum_remap_name = self.inputs_cfg.get('INPUTS', 'covden_sum_remap')
        except:
            self.covden_sum_remap_name = None
            
        try:
            self.covden_win_remap_name = self.inputs_cfg.get('INPUTS', 'covden_win_remap')
        except:
            self.covden_win_remap_name = None
            
        try:
            self.snow_intcp_remap_name = self.inputs_cfg.get('INPUTS', 'snow_intcp_remap')
        except:
            self.snow_intcp_remap_name = None
            
        try:
            self.srain_intcp_remap_name = self.inputs_cfg.get('INPUTS', 'srain_intcp_remap')
        except:
            self.srain_intcp_remap_name = None
            
        try:
            self.wrain_intcp_remap_name = self.inputs_cfg.get('INPUTS', 'wrain_intcp_remap')
        except:
            self.wrain_intcp_remap_name = None
            
        try: 
            self.root_depth_remap_name = self.inputs_cfg.get('INPUTS', 'root_depth_remap')
        except:
            self.root_depth_remap_name = None
        
        # Check that remap folder is valid
        if not os.path.isdir(self.remap_ws):
            logging.error('\nERROR: Remap folder does not exist\n')
            sys.exit()
                
        
    def read_DEM_parameters(self):
        """
        Read the DEM properties from the config file
        """
        
        self.dem_orig_path = self.inputs_cfg.get('INPUTS', 'dem_orig_path')
        
        # Resampling method 'BILINEAR', 'CUBIC', 'NEAREST'
        self.dem_proj_method = self.inputs_cfg.get('INPUTS', 'dem_projection_method').upper()
        self.dem_cs = self.inputs_cfg.getint('INPUTS', 'dem_cellsize')
    
        #
        self.reset_dem_adj_flag = self.inputs_cfg.getboolean('INPUTS', 'reset_dem_adj_flag')
        self.dem_adj_copy_field = self.inputs_cfg.get('INPUTS', 'dem_adj_copy_field')
    
        # Use PRISM temperature to set Jensen-Haise coefficient
        self.calc_prism_jh_coef_flag = self.inputs_cfg.getboolean(
            'INPUTS', 'calc_prism_jh_coef_flag')
    
        # Calculate flow accumulation weighted elevation
        self.calc_flow_acc_dem_flag = self.inputs_cfg.getboolean(
            'INPUTS', 'calc_flow_acc_dem_flag')
        
        # Check that either the original DEM raster exists
        if not arcpy.Exists(self.dem_orig_path):
            logging.error(
                '\nERROR: DEM ({0}) raster does not exist\n'.format(self.dem_orig_path))
            sys.exit()
        
        # Check other inputs
        if self.dem_cs <= 0:
            logging.error('\nERROR: DEM cellsize must be greater than 0')
            sys.exit()
            
        dem_proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if self.dem_proj_method not in dem_proj_method_list:
            logging.error('\nERROR: DEM projection method must be: {0}'.format(
                ', '.join(dem_proj_method_list)))
            sys.exit()
            
        if self.reset_dem_adj_flag:
            logging.warning('\nWARNING: All values in {0} will be overwritten'.format(
                self.dem_adj_field))
            raw_input('  Press ENTER to continue')
            
        # Check that dem_adj_copy_field exists
        if len(arcpy.ListFields(self.polygon_path, self.dem_adj_copy_field)) == 0:
            logging.error('\nERROR: dem_adj_copy_field {0} does not exist\n'.format(
                self.dem_adj_copy_field))
            sys.exit()
        
        
    def read_veg_parameters(self):
        """
        Read the vegetation parameters from the config file
        """
        
        # Landfire Vegetation Type
        self.veg_type_orig_path = self.inputs_cfg.get('INPUTS', 'veg_type_orig_path')
        self.veg_type_cs = self.inputs_cfg.getint('INPUTS', 'veg_type_cellsize')
        try:
            self.veg_type_field = self.inputs_cfg.get('INPUTS', 'veg_type_field')
        except:
            self.veg_type_field = None
    
        # Landfire Vegetation Cover
        self.veg_cover_orig_path = self.inputs_cfg.get('INPUTS', 'veg_cover_orig_path')
        self.veg_cover_cs = self.inputs_cfg.getint('INPUTS', 'veg_cover_cellsize')
        
        # Check that either the original vegetation raster exist
        if not arcpy.Exists(self.veg_cover_orig_path):
            logging.error(
                '\nERROR: Vegetation cover raster does not exist')
            sys.exit()
        if not arcpy.Exists(self.veg_type_orig_path):
            logging.error(
                '\nERROR: Vegetation type raster does not exist')
            sys.exit()
            
        # Vegetation cover can be set from another field in the raster
        # This is mostly for US_120EVT
        if not self.veg_type_field:
            logging.info('\n  Using VALUE field to set vegetation type')
            veg_type_field = 'VALUE'
        elif len(arcpy.ListFields(self.veg_type_orig_path, self.veg_type_field)) == 0:
            logging.info(
                ('  veg_type_field {0} does not exist\n  Using VALUE ' +
                 'field to set vegetation type').format(veg_type_field))
            veg_type_field = 'VALUE'
        elif arcpy.ListFields(self.veg_type_orig_path, self.veg_type_field)[0].type not in ['Integer', 'SmallInteger']:
            logging.info(
                ('  veg_type_field {0} is not an integer type\n  Using VALUE ' +
                 'field to set vegetation type').format(self.veg_type_field))
            self.veg_type_field = 'VALUE'
            
        # Check other inputs
        if self.veg_type_cs <= 0:
            logging.error('\nERROR: Veg. type cellsize must be greater than 0')
            sys.exit()
        if self.veg_cover_cs <= 0:
            logging.error('\nERROR: Veg. cover cellsize must be greater than 0')
            sys.exit()
            
    def read_soil_parameters(self):
        """
        Read the soil parameters from the config file
        """
        
        self.soil_orig_ws  = self.inputs_cfg.get('INPUTS', 'soil_orig_folder')
        self.awc_name      = self.inputs_cfg.get('INPUTS', 'awc_name')
        self.clay_pct_name = self.inputs_cfg.get('INPUTS', 'clay_pct_name')
        self.sand_pct_name = self.inputs_cfg.get('INPUTS', 'sand_pct_name')
        # silt_pct_name = self.inputs_cfg.get('INPUTS', 'silt_pct_name')
        self.soil_proj_method = 'NEAREST'
        self.soil_cs = self.inputs_cfg.getint('INPUTS', 'soil_cellsize')
        self.fill_soil_nodata_flag = self.inputs_cfg.getboolean(
            'INPUTS', 'fill_soil_nodata_flag')
        
        self.soil_pct_flag = self.inputs_cfg.getboolean('INPUTS', 'soil_pct_flag')
        self.moist_init_ratio = self.inputs_cfg.getfloat('INPUTS', 'moist_init_ratio')
        self.rechr_init_ratio = self.inputs_cfg.getfloat('INPUTS', 'rechr_init_ratio')
    
        # Use Ksat to calculate ssr2gw_rate and slowcoef_lin
        # calc_ssr2gw_rate_flag = hru.inputs_cfg.getboolean(
        #    'INPUTS', 'calc_ssr2gw_rate_flag')
        # calc_slowcoef_flag = hru.inputs_cfg.getboolean(
        #    'INPUTS', 'calc_slowcoef_flag')
        # if calc_ssr2gw_rate_flag or calc_slowcoef_flag:
        self.ksat_name = self.inputs_cfg.get('INPUTS', 'ksat_name')
    
        # Clip root depth to soil depth
        self.clip_root_depth_flag = self.inputs_cfg.getboolean(
            'INPUTS', 'clip_root_depth_flag')
        if self.clip_root_depth_flag:
            self.soil_depth_name = self.inputs_cfg.get('INPUTS', 'soil_depth_name')
        
        # All of the soil rasters must exist
        self.awc_orig_path = os.path.join(self.soil_orig_ws, self.awc_name)
        self.clay_pct_orig_path = os.path.join(self.soil_orig_ws, self.clay_pct_name)
        self.sand_pct_orig_path = os.path.join(self.soil_orig_ws, self.sand_pct_name)
        # silt_orig_path = os.path.join(soil_orig_ws, silt_pct_name)
        
        # if calc_ssr2gw_rate_flag or calc_slowcoef_flag:
        self.ksat_orig_path = os.path.join(self.soil_orig_ws, self.ksat_name)
        if self.clip_root_depth_flag:
            self.soil_depth_path = os.path.join(self.soil_orig_ws, self.soil_depth_name)
    
        # Check soil init ratios
        if self.moist_init_ratio < 0 or self.moist_init_ratio > 1:
            logging.error('\nERROR: Soil moist_init_ratio must be between 0 & 1')
            sys.exit()
        if self.rechr_init_ratio < 0 or self.rechr_init_ratio > 1:
            logging.error('\nERROR: Soil rechr_init_ratio must be between 0 & 1')
            sys.exit()
    
    
        # Check that either the original or projected/clipped raster exists
        if not arcpy.Exists(self.awc_orig_path):
            logging.error('\nERROR: AWC raster does not exist')
            sys.exit()
        if not arcpy.Exists(self.clay_pct_orig_path):
            logging.error('\nERROR: Clay raster does not exist')
            sys.exit()
        if not arcpy.Exists(self.sand_pct_orig_path):
            logging.error('\nERROR: Sand raster does not exist')
            sys.exit()
        # if not arcpy.Exists(silt_orig_path):
        #    logging.error('\nERROR: Silt raster does not exist')
        #    sys.exit()
        # if ((calc_ssr2gw_rate_flag or calc_slowcoef_flag) and
        #    not arcpy.Exists(ksat_orig_path)):
        if not arcpy.Exists(self.ksat_orig_path):
            logging.error('\nERROR: Ksat raster does not exist')
            sys.exit()
        if self.clip_root_depth_flag and not arcpy.Exists(self.soil_depth_orig_path):
            logging.error('\nERROR: Soil depth raster does not exist')
            sys.exit()
    
        # Check other inputs
        if self.soil_cs <= 0:
            logging.error('\nERROR: soil cellsize must be greater than 0')
            sys.exit()
        soil_proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if self.soil_proj_method.upper() not in soil_proj_method_list:
            logging.error('\nERROR: Soil projection method must be: {0}'.format(
                ', '.join(soil_proj_method_list)))
            sys.exit()
        
        
    def read_impervious_parameters(self):
        """
        Read the impervious parameters from the config file
        """
        
        #
        self.imperv_orig_path = self.inputs_cfg.get('INPUTS', 'impervious_orig_path')
        # imperv_proj_method = inputs_cfg.get('INPUTS', 'impervious_projection_method')
        self.imperv_proj_method = 'NEAREST'
        self.imperv_cs = self.inputs_cfg.getint('INPUTS', 'impervious_cellsize')
        self.imperv_pct_flag = self.inputs_cfg.getboolean('INPUTS', 'impervious_pct_flag')
        
        # Impervious raster must exist
        if not arcpy.Exists(self.imperv_orig_path):
            logging.error('\nERROR: Impervious raster does not exist')
            sys.exit()
    
        # Check other inputs
        if self.imperv_cs <= 0:
            logging.error('\nERROR: soil cellsize must be greater than 0')
            sys.exit()
        imperv_proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if self.imperv_proj_method.upper() not in imperv_proj_method_list:
            logging.error('\nERROR: Impervious projection method must be: {0}'.format(
                ', '.join(imperv_proj_method_list)))
            sys.exit()

    def read_prism_parameters(self):
        """
        Read the PRISM parameters from the config file
        """
        
        # PRISM
        self.prism_ws = self.inputs_cfg.get('INPUTS', 'prism_folder')
        self.prism_proj_method = self.inputs_cfg.get('INPUTS', 'prism_projection_method')
        self.prism_cs = self.inputs_cfg.getint('INPUTS', 'prism_cellsize')
        self.calc_prism_jh_coef_flag = self.inputs_cfg.getboolean(
            'INPUTS', 'calc_prism_jh_coef_flag')
            
        # Check that PRISM folder is valid
        if not os.path.isdir(self.prism_ws):
            logging.error(
                '\nERROR: PRISM folder ({0}) does not exist'.format(self.prism_ws))
            sys.exit()
            
        proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if self.prism_proj_method.upper() not in proj_method_list:
            logging.error('\nERROR: PRISM projection method must be: {0}'.format(
                ', '.join(proj_method_list)))
            sys.exit()
        logging.debug('  Projection method:    {0}'.format(
            self.prism_proj_method.upper()))
    
        # Check other inputs
        if self.prism_cs <= 0:
            logging.error('\nERROR: PRISM cellsize must be greater than 0\n')
            sys.exit()
            
    def read_stream_parameters(self):
        """
        Read the stream parameters from the config file
        """
        
        self.streams_path = self.inputs_cfg.get('INPUTS', 'streams_path')
        if not os.path.isfile(self.streams_path):
            logging.error(
                ('\nERROR: Stream shapefiles does not exist' +
                 '\nERROR:   {0}').format(
                     self.streams_path))
            sys.exit()


def next_row_col(flow_dir, cell):
    """"""
    i_next, j_next = cell
    # Upper left cell is 0,0
    if flow_dir in [1, 2, 128]:
        i_next += 1
    elif flow_dir in [8, 16, 32]:
        i_next -= 1
    if flow_dir in [2, 4, 8]:
        j_next += 1
    elif flow_dir in [32, 64, 128]:
        j_next -= 1
    return i_next, j_next


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


def add_field_func(hru_param_path, field_name, field_type='DOUBLE'):
    """"""
    while not arcpy.ListFields(hru_param_path, field_name):
        logging.info('  Field: {0}'.format(field_name))
#         try:
        arcpy.AddField_management(hru_param_path, field_name, field_type)
#         except:
#             pass
#         sleep(0.5)


def transform_func(spat_ref_a, spat_ref_b):
    """"""
    # Set preferred transforms
    if ((spat_ref_a.GCS.name == 'GCS_WGS_1984' and
         spat_ref_b.GCS.name == 'GCS_North_American_1983') or
        (spat_ref_a.GCS.name == 'GCS_North_American_1983' and
         spat_ref_b.GCS.name == 'GCS_WGS_1984')):
        return 'NAD_1983_To_WGS_1984_5'
    else:
        return None
        # return '#'


def valid_raster_func(raster_path, raster_name, hru_param, cs=10):
    """This will check for matching spat. ref., snap point, and cellsize

    Assume raster is not valid unless it passes all checks
    """
    if not arcpy.Exists(raster_path):
        return False
    logging.debug('\nReading existing {0} raster'.format(raster_name))
    raster_obj = Raster(raster_path)
    raster_sr = raster_obj.spatialReference
    raster_extent = raster_obj.extent
    raster_cs = raster_obj.meanCellWidth
    logging.debug('  {0} spat. ref.: {1}'.format(
        raster_name, raster_sr.name))
    logging.debug('  {0} GCS:        {1}'.format(
        raster_name, raster_sr.GCS.name))
    logging.debug('  {0} extent:     {1}'.format(
        raster_name, extent_string(raster_extent)))
    logging.debug('  {0} cellsize:   {1}'.format(raster_name, raster_cs))
    ref_pnt = arcpy.Point(hru_param.ref_x, hru_param.ref_y)
    if raster_sr.name != hru_param.sr.name:
        logging.error(
            ('\nERROR: The {0} spatial reference does not match ' +
             'the hru_param spatial reference').format(raster_name))
        return False
    elif not snapped(raster_extent, ref_pnt, hru_param.cs):
        logging.error(
            '\nWARNING: The {0} is not snapped to the hru_param'.format(
                raster_name))
        return False
    elif raster_cs != cs:
        logging.error(
            '\nERROR: The {0} needs to have a {1}m cellsize'.format(
                raster_name, cs))
        return False
    elif not raster_extent.contains(hru_param.extent):
        logging.error(
            '\nERROR: The {0} extent is too small'.format(raster_name))
        return False
    else:
        return True


def zonal_stats_func(zs_dict, polygon_path, hru_param,
                     nodata_value=-999, default_value=0):
    """
    Calculate zonal statistics for each HRU
    
    - Get subset of HRU polygons
    - Zonal stats by table based on polygons
    - Add the zonal stats back to the HRU polygons
    """
    
    for zs_field, (raster_path, zs_stat) in sorted(zs_dict.items()):
        logging.info('  {0}: {1}'.format(zs_field, zs_stat))
        logging.info('    {0}'.format(raster_path))
        # Check inputs
        zs_stat_list = ['MEAN', 'MINIMUM', 'MAXIMUM', 'MEDIAN', 'MAJORITY', 'SUM']
        zs_field_list = arcpy.ListFields(polygon_path, zs_field)
        if zs_stat not in zs_stat_list:
            sys.exit()
        elif len(zs_field_list) == 0:
            logging.error(
                '\nERROR: Zonal stats field {0} doesn\'t exist'.format(zs_field))
            sys.exit()

    # Check that the shapefiles have a spatial reference
    if arcpy.Describe(polygon_path).spatialReference.name == 'Unknown':
        logging.error(
            '\nERROR: HRU centroids  is not projected (i.e. does not have a prj file)')
        sys.exit()
        
#     if arcpy.Describe(point_path).spatialReference.name == 'Unknown':
#         logging.error(
#             '\nERROR: HRU centroids does not appear to be projected (or does not have a prj file)' +
#             '\nERROR: Try deleting the centroids (i.e. "_label.shp") and ' +
#             'rerunning hru_parameters.py\n')
#         sys.exit()

    # Check that ORIG_FID is in point_path (HRU centroids)
    if len(arcpy.ListFields(polygon_path, hru_param.fid_field)) == 0:
        logging.error(
            ('\nERROR: HRU polygons does not have the field: {0}' +
             '\nERROR: Try deleting the centroids (i.e. "_label.shp") and ' +
             'rerunning hru_parameters.py\n').format(hru_param.fid_field))
        sys.exit()

    # Check for duplicate ORIG_FID values
    hru_param_count = int(arcpy.GetCount_management(polygon_path).getOutput(0))
    if field_duplicate_check(polygon_path, hru_param.fid_field, hru_param_count):
        logging.error(
            ('\nERROR: There are duplicate {0} values\n').format(hru_param.fid_field))
        sys.exit()
        
    # DEADBEEF - remove once field_duplicate_check() is full developed
    # fid_list = [r[0] for r in arcpy.da.SearchCursor(point_path, [hru_param.fid_field])]
    # if len(fid_list) != len(set(fid_list)):
    #    logging.error(
    #        ('\nERROR: There are duplicate {0} values\n').format(hru_param.fid_field))
    #    sys.exit()

    # Create memory objects
    polygon_subset_path = os.path.join('in_memory', 'polygon_subset')
    hru_raster_path = os.path.join('in_memory', 'hru_raster')
#     polygon_subset_path = os.path.join(os.path.join(hru_param.param_ws, 'dem_rasters'), 'polygon_subset.shp')
#     hru_raster_path = os.path.join(os.path.join(hru_param.param_ws, 'dem_rasters'), 'hru_raster')
    # point_subset_path = os.path.join(env.scratchWorkspace, 'point_subset.shp')
    # hru_raster_path = os.path.join(env.scratchWorkspace, 'hru_raster.img')
    # Set environment parameters for polygon to raster conversion
    env.extent = hru_param.extent
    env.outputCoordinateSystem = polygon_path
    # env.cellSize = hru_param.cs

    # Only ~65536 objects can be processed by zonal stats
    block_size = 65000
    for i, x in enumerate(xrange(0, hru_param_count, block_size)):
        logging.info('  FIDS: {0}-{1}'.format(x, x + block_size))
        
        # Select a subset of the cell centroids
        logging.debug('    Selecting FID subset')
        subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
            hru_param.fid_field, x, x + block_size)
        arcpy.Select_analysis(
            polygon_path, polygon_subset_path, subset_str)
        
        # Convert points subset to raster
#         logging.debug('    Converting shapefile to raster')
#         arcpy.FeatureToRaster_conversion(
#             polygon_subset_path, hru_param.fid_field,
#             hru_raster_path, hru_param.dem_cs)

        # Zonal stats
        logging.debug('    Calculating zonal stats')
        data_dict = defaultdict(dict)
        for zs_field, (raster_path, zs_stat) in sorted(zs_dict.items()):
            zs_name = '{0}_{1}'.format(zs_field.upper(), i)
            logging.info('    {0}: {1}'.format(zs_stat.upper(), zs_name))
            # For some reason with 10.2, ZS doesn't work with cs at HRU cs
            env.cellSize = Raster(raster_path).meanCellWidth
            
            # Calculate zonal statistics
            zs_table = os.path.join('in_memory', zs_name)
#             zs_table = os.path.join(hru_param.param_ws, 'dem_rasters', zs_name)
            zs_obj = ZonalStatisticsAsTable(
                polygon_subset_path, hru_param.fid_field, raster_path,
                zs_table, 'DATA', zs_stat.upper())

            # Read values from points
            logging.debug('    Reading values from zs table')
            # Fields 1 & 4 are the 'Value' (ORIG_FID) and the stat (SUM, MEAN, etc)
            fields = [
                f.name for f_i, f in enumerate(arcpy.ListFields(zs_table))
                if f_i in [1, 4]]
            for row in arcpy.da.SearchCursor(zs_table, fields):
                # Set NoData value for cells that are entirely NoData
                if row[1] is None:
                    data_dict[int(row[0])][zs_field] = nodata_value
                else:
                    data_dict[int(row[0])][zs_field] = float(row[1])
            arcpy.Delete_management(zs_obj)
            arcpy.Delete_management(zs_table)
            del zs_table, zs_obj, fields

        # Write values to polygon
        logging.info('    Writing values to polygons')
        zs_fields = sorted(zs_dict.keys())
        fields = zs_fields + [hru_param.fid_field]
        with arcpy.da.UpdateCursor(polygon_path, fields, subset_str) as u_cursor:
            for row in u_cursor:
                # Create an empty dictionary if FID does not exist
                # Missing FIDs did not have zonal stats calculated
                row_dict = data_dict.get(int(row[-1]), None)
                for i, zs_field in enumerate(zs_fields):
                    # If stats were calculated for only some parameters,
                    #   then set missing parameter value to nodata value (-999)
                    if row_dict:
                        try:
                            row[i] = row_dict[zs_field]
                        except KeyError:
                            row[i] = nodata_value
                    # Otherwise, if no stats were calculated,
                    #   reset value to 0 (shapefile default)
                    else:
                        row[i] = default_value
                u_cursor.updateRow(row)

        # Cleanup
        del data_dict
        if arcpy.Exists(polygon_subset_path):
            arcpy.Delete_management(polygon_subset_path)
        if arcpy.Exists(hru_raster_path):
            arcpy.Delete_management(hru_raster_path)

    arcpy.ClearEnvironment('extent')
    arcpy.ClearEnvironment('outputCoordinateSystem')
    arcpy.ClearEnvironment('cellSize')


def field_duplicate_check(table_path, field_name, n=None):
    """Check if there are duplicate values in a shapefile field

    For now assume table_path is actually a shapefile that can be read
        with arcpy.da.SearchCursor()

    Args:
        table_path (str): File path of the table to search
        field_name (str): Field/column name to search

    Returns:
        bool: True if there are duplicate values in the field, False otherwise
    """

    # Eventually check that field is in table
    field_obj = arcpy.ListFields(table_path, field_name)[0]

    if n is None:
        n = int(arcpy.GetCount_management(table_path).getOutput(0))
    n32_max = 3000000
    logging.debug('\n  Testing for duplicate values')
    logging.debug('    field:    {}'.format(field_name))
    logging.debug('    features: {}'.format(n))
    logging.debug('    n32_max:  {}'.format(n32_max))
    logging.debug('    2**32:    {}'.format(2**32))
    logging.debug('    maxsize:  {}'.format(sys.maxsize))

    if sys.maxsize > 2**32 or n < n32_max:
        logging.debug('    Reading values')
        # If 64-bit or row count is low, read all values into memory
        fid_list = [r[0] for r in arcpy.da.SearchCursor(table_path, [field_name])]
        if len(fid_list) != len(set(fid_list)):
            return True
        else:
            logging.debug('    No duplicates')
            return False
    elif field_obj.type in ['Integer', 'SmallInteger']:
        # This approach will only work with integers
        block_size = 500000
        fid_ranges = []
        for i, x in enumerate(xrange(0, n, block_size)):
            logging.debug('    FIDS: {0}-{1}'.format(x, x + block_size))
            subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
                arcpy.Describe(table_path).OIDFieldName, x, x + block_size)

            # Don't sort here since it gets sorted in group_ranges()
            fid_list = [r[0] for r in arcpy.da.SearchCursor(
                table_path, [field_name], subset_str)]

            # Return True if there are duplicates in a subset
            if len(fid_list) != len(set(fid_list)):
                return True

            # Group consecutive values into ranges
            fid_new_ranges = list(group_ranges(fid_list))

            # Check if any of the new ranges overlap each other
            # Skip for now since subset duplicates were checked above
            # if ranges_overlap(fid_new_ranges):
            #    return True

            # Check if any of the new ranges overlap the existing ranges
            if ranges_overlap(fid_ranges + fid_new_ranges):
                return True

            # Merge subset ranges into main range list
            fid_ranges = list(merge_ranges(fid_ranges + fid_new_ranges))
            del fid_new_ranges

        logging.debug('    FID ranges: {}'.format(fid_ranges))
        logging.debug('    No duplicates')
        return False
    else:
        # For now, assume there are not duplicates if the values can't
        #   all be read in
        logging.debug(
            '    Assuming no duplicates since field type is not integer, \n' +
            '      table contains more than {} rows, and 32-bit Python is \n' +
            '      limited to 2 GB of memory.')
        return False
    # DEADBEEF - I'm not sure this will actually work on a really large table
    # else:
    #    duplicate_flag = False
    #    fid_prev = None
    #    cursor = arcpy.SearchCursor(
    #        table_path, fields=field_name, sort_fields=field_name+' A')
    #    for row in cursor:
    #        fid = row.getValue(field_name)
    #        if fid == fid_prev and fid_prev is not None:
    #            duplicate_flag = True
    #            break
    #        else:
    #            fid_prev = fid
    #    del cursor, row
    #    logging.debug('    No duplicates')
    #    return duplicate_flag


def group_ranges(input_list):
    """Group

    Copied from:
    http://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list

    Args:
        input_list (list): list of numbers to group into ranges

    Yields
         tuple: pairs (min, max)
    """
    for k, g in itertools.groupby(enumerate(sorted(input_list)), lambda (i, x):i-x):
        group = map(itemgetter(1), g)
        yield group[0], group[-1]


def merge_ranges(ranges):
    """Merge overlapping and adjacent integer ranges

    Yield the merged ranges in order
    The argument must be an iterable of pairs (start, stop).

    Copied from:
    http://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []

    Args:
        ranges (list): Iterable of pairs (min, max)

    Yields:
        tuple: pairs (min, max)
    """
    ranges = iter(sorted(ranges))
    current_start, current_stop = next(ranges)
    for start, stop in ranges:
        if start > (current_stop + 1):
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop


def ranges_overlap(ranges):
    """Test if ranges overlap

    Args:
        ranges (list): Iterable of pairs (min, max)
    Returns:
         bool: True if ranges overlap each other, False otherwise
    """
    for r1, r2 in itertools.combinations(ranges, 2):
        if r1[1] > r2[0] and r1[0] < r2[1]:
            return True
    return False


def extent_string(extent_obj):
    """"""
    return ' '.join(str(extent_obj).split()[:4])
    # return ' '.join(['{0:.4f}'.format(s) for s in str(extent_obj).split()[:4]])


def round_extent(extent_obj, n=10):
    """"""
    return arcpy.Extent(
        round(extent_obj.XMin, n), round(extent_obj.YMin, n),
        round(extent_obj.XMax, n), round(extent_obj.YMax, n))


# This adjusts one extent to a snap point
# This is similar to the GDAL implementation
def adjust_extent_to_snap(extent_obj, snap_pnt, cs, method='EXPAND',
                          integer_flag=True):
    """"""
    if method.upper() == 'ROUND':
        extent_xmin = math.floor(
            (extent_obj.XMin - snap_pnt.X) / cs + 0.5) * cs + snap_pnt.X
        extent_ymin = math.floor(
            (extent_obj.YMin - snap_pnt.Y) / cs + 0.5) * cs + snap_pnt.Y
        extent_xmax = math.floor(
            (extent_obj.XMax - snap_pnt.X) / cs + 0.5) * cs + snap_pnt.X
        extent_ymax = math.floor(
            (extent_obj.YMax - snap_pnt.Y) / cs + 0.5) * cs + snap_pnt.Y
    elif method.upper() == 'EXPAND':
        extent_xmin = math.floor(
            (extent_obj.XMin - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymin = math.floor(
            (extent_obj.YMin - snap_pnt.Y) / cs) * cs + snap_pnt.Y
        extent_xmax = math.ceil(
            (extent_obj.XMax - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymax = math.ceil(
            (extent_obj.YMax - snap_pnt.Y) / cs) * cs + snap_pnt.Y
    elif method.upper() == 'SHRINK':
        extent_xmin = math.ceil(
            (extent_obj.XMin - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymin = math.ceil(
            (extent_obj.YMin - snap_pnt.Y) / cs) * cs + snap_pnt.Y
        extent_xmax = math.floor(
            (extent_obj.XMax - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymax = math.floor(
            (extent_obj.YMax - snap_pnt.Y) / cs) * cs + snap_pnt.Y
    if integer_flag:
        return arcpy.Extent(
            int(round(extent_xmin, 0)), int(round(extent_ymin, 0)),
            int(round(extent_xmax, 0)), int(round(extent_ymax, 0)))
    else:
        return arcpy.Extent(extent_xmin, extent_ymin, extent_xmax, extent_ymax)


def buffer_extent_func(extent_obj, extent_buffer):
    """"""
    # extent_obj = arcpy.Describe(extent_feature).extent
    extent_xmin = extent_obj.XMin-extent_buffer
    extent_ymin = extent_obj.YMin-extent_buffer
    extent_xmax = extent_obj.XMax+extent_buffer
    extent_ymax = extent_obj.YMax+extent_buffer
    return arcpy.Extent(
        extent_xmin, extent_ymin, extent_xmax, extent_ymax)


# Check if rasters are aligned to snap_raster
# Check if rasters have same cellsize as snap_raster
def snapped(extent_obj, snap_pnt, cs):
    """"""
    if (((snap_pnt.X - extent_obj.XMin) % cs == 0) and
        ((snap_pnt.X - extent_obj.XMax) % cs == 0) and
        ((snap_pnt.Y - extent_obj.YMin) % cs == 0) and
        ((snap_pnt.Y - extent_obj.YMax) % cs == 0)):
        return True
    else:
        return False


def get_ini_file(workspace, ini_re, function_str='function'):
    """"""
    # Get ini file name
    ini_file_list = build_file_list(workspace, ini_re)
    # Filter field list ini file
    ini_file_list = [
        item for item in ini_file_list if 'field_list.ini' not in item]
    if len(ini_file_list) == 1:
        config_filepath = ini_file_list[0]
    elif len(ini_file_list) > 1:
        ini_file_len_max = max([len(item) for item in ini_file_list])
        # An ini file was not passed as an arguement to the script
        # Look for ini files in the working directory
        # If only one, use it
        # If more, let user pick from list)
        # If none, error out
        print('\nThere is more than one INI file present in the folder')
        print('  {0:2s}  {1}'.format('# ', 'INI File'))
        print('  {0:2s}  {1}'.format('==', '='*ini_file_len_max))
        for i, ini_file in enumerate(ini_file_list):
            print('  {0:2d}  {1}'.format(i, ini_file))
        config_filepath = None
        while not config_filepath:
            usr_input = raw_input('\nPlease select an INI file to use: ')
            try:
                ini_file_index = int(usr_input)
                config_filepath = ini_file_list[ini_file_index]
            except (ValueError, IndexError):
                pass
        print('  Using {0}\n'.format(config_filepath))
        del ini_file_len_max, usr_input
    else:
        print('\nERROR: No suitable ini files were found')
        print('ERROR: Please set input file when calling {0}'.format(function_str))
        print('ERROR: For example: test.py test.ini\n')
        sys.exit()
    config_filename = os.path.basename(config_filepath)
    print('{0:<20s} {1}'.format('INI File Name:', config_filename))
    return config_filepath


def get_param(param_str, param_default, config, section='INPUTS'):
    """"""
    param_type = type(param_default)
    try:
        if param_type is float:
            param_value = config.getfloat('INPUTS', param_str)
        elif param_type is int:
            param_value = config.getint('INPUTS', param_str)
        elif param_type is bool:
            param_value = config.getboolean('INPUTS', param_str)
        elif param_type is list or param_type is tuple:
            param_value = [
                # i for i in re.split('\W+', config.get('INPUTS', param_str)) if i]
                i.strip() for i in config.get('INPUTS', param_str).split(',')
                if i.strip()]
        elif param_type is str or param_default is None:
            param_value = config.get('INPUTS', param_str)
            if param_value.upper() == 'NONE':
                param_value = None
        else:
            logging.error('ERROR: Unknown Input Type: {0}'.format(param_type))
            sys.exit()
    except:
        param_value = param_default
        if param_type is str and param_value.upper() == 'NONE':
            param_value = None
        logging.warning('  NOTE: {0} = {1}'.format(param_str, param_value))
    return param_value


def build_file_list(ws, test_re, test_other_re=None):
    """"""
    if test_other_re is None:
        test_other_re = re.compile('a^')
    if os.path.isdir(ws):
        return sorted([os.path.join(ws, item) for item in os.listdir(ws)
                       if (os.path.isfile(os.path.join(ws, item)) and
                           test_re.match(item) or test_other_re.match(item))])
    else:
        return []


def get_prism_data_name():
    """"""
    #  Get PRISM data name
    data_name_dict = dict()
    data_name_dict[1] = 'PPT'
    data_name_dict[2] = 'TMAX'
    data_name_dict[3] = 'TMIN'
    data_name_dict[4] = 'ALL'
    print '\nPlease select which PRISM product(s) to calculate'
    print '  {0:2s}  {1}'.format('# ', 'PRISM Data')
    print '  {0:2s}  {1}'.format('==', '==========')
    for i, data_name in sorted(data_name_dict.items()):
        print '  {0:2d}  {1}'.format(i, data_name)
    data_name = None
    while not data_name:
        usr_input = raw_input(
            '\nPlease select a PRISM data product to calculate: ')
        try:
            data_name_index = int(usr_input)
            data_name = data_name_dict[data_name_index]
        except (ValueError, IndexError):
            pass
    print '  Using {0}\n'.format(data_name)
    return data_name


def project_hru_extent_func(hru_extent, hru_sr,
                            target_extent, target_cs, target_sr):
    """"""
    
    # uses an abitrary cell size of 300 to create a polygon of the extent,
    # one could also set step to ensure that there are a certain number of steps
    # along the extent boundary
    hru_cs=300
    
    logging.debug('  Projecting extent')
    logging.debug('  HRU Extent:   {0}'.format(extent_string(hru_extent)))
#     logging.debug('  HRU cellsize: {0}'.format(hru_cs))
    logging.debug('  HRU spatref:  {0}'.format(hru_sr.name))
    logging.debug('  Target snap:     {0}'.format(target_extent.lowerLeft))
    logging.debug('  Target cellsize: {0}'.format(target_cs))
    logging.debug('  Target spatref:  {0}'.format(target_sr.name))
    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    # Project the HRU extent to the raster spatial reference
    hru_corners = [
        [hru_extent.XMin, hru_extent.YMax],
        [hru_extent.XMax, hru_extent.YMax],
        [hru_extent.XMax, hru_extent.YMin],
        [hru_extent.XMin, hru_extent.YMin],
        [hru_extent.XMin, hru_extent.YMax]]
    
    # Add points between corners
    hru_points = []
    for point_a, point_b in zip(hru_corners[:-1], hru_corners[1:]):
        steps = float(max(
            abs(point_b[0] - point_a[0]),
            abs(point_b[1] - point_a[1]))) / hru_cs
        for x, y in zip(np.linspace(point_a[0], point_b[0], steps + 1),
                        np.linspace(point_a[1], point_b[1], steps + 1)):
            hru_points.append(arcpy.Point(x,y))
            
    # Project all points to output spatial reference and get projected extent
    transform = transform_func(hru_sr, target_sr)
    if transform:
        projected_extent = arcpy.Polygon(
            arcpy.Array(hru_points), hru_sr).projectAs(
                target_sr, transform).extent
    else:
        projected_extent = arcpy.Polygon(
            arcpy.Array(hru_points), hru_sr).projectAs(target_sr).extent
    logging.debug('  Projected Extent: {0}'.format(
        extent_string(projected_extent)))
    
    # Adjust extent to match snap
#     projected_extent = adjust_extent_to_snap(
#         projected_extent, target_extent.lowerLeft, target_cs, 'EXPAND', False)
#     logging.debug('  Snapped Extent:   {0}'.format(
#         extent_string(projected_extent)))
    
    # Buffer extent 4 input cells
    # projected_extent = buffer_extent_func(projected_extent, 4 * target_cs)
#     projected_extent = buffer_extent_func(
#         projected_extent, 4 * max(target_cs, hru_cs))
#     logging.debug('  Buffered Extent:  {0}'.format(
#         extent_string(projected_extent)))
    return projected_extent


def project_raster_func(input_raster, output_raster, output_sr,
                        proj_method, input_cs, transform_str,
                        input_sr, hru_param):
    """"""
    # Input raster can be a raster object or a raster path
    # print isinstance(input_raster, Raster), isinstance(input_raster, str)
    try:
        input_extent = Raster(input_raster).extent
    except:
        input_extent = input_raster.extent
        
    # DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    # Clip the input raster with the projected HRU extent first
    # Project extent from "output" to "input" to get clipping extent
    proj_extent = project_hru_extent_func(
        hru_param.extent, output_sr,
        input_extent, input_cs, input_sr)
    # clip_path = output_raster.replace('.img', '_clip.img')
    clip_path = os.path.join('in_memory', 'clip_raster')
#     clip_path = os.path.join(hru_param.param_ws, 'clip_raster')
    env.extent = proj_extent
    arcpy.Clip_management(
        input_raster, ' '.join(str(proj_extent).split()[:4]), clip_path)
    arcpy.ClearEnvironment('extent')
    
    # Then project the clipped raster
#     arcpy.ProjectRaster_management(
#         clip_path, output_raster, output_sr, proj_method.upper(), input_cs,
#         transform_str, reg_point, input_sr)
    arcpy.ProjectRaster_management(in_raster=clip_path, out_raster=output_raster, 
                                   out_coor_system=output_sr, resampling_type=proj_method.upper(),
                                   cell_size=input_cs, geographic_transform=transform_str, 
                                   in_coor_system=input_sr)
    
    # Cleanup
    arcpy.Delete_management(clip_path)


def cell_area_func(hru_param_path, area_field):
    """"""
    arcpy.CalculateField_management(
        hru_param_path, area_field, '!SHAPE.AREA@acres!', 'PYTHON')


def zone_by_area_func(zone_path, zone_field, zone_value, hru_param_path,
                      hru_param, hru_area_field='HRU_AREA',
                      zone_area_field=None, area_pct=50):
    """Flag cells that are inside a feature based on an area weighting

    Set values that are in zone, but don't reset values that are out of zone

    Args:
        zone_path (str):
        zone_field (str):
        zone_value (int):
        hru_param_path (str):
        hru_param: class:`support_functions.HRUParameters`
        hru_area_field (str):
        zone_area_field (str):
        area_pct ():

    Returns:
        None
    """
    zone_value_field = 'ZONE_VALUE'
    int_area_field = 'INT_AREA'
    int_pct_field = 'INT_PCT'
    # Need to set zone value into a field before intersect
    # If zone_value is FID, add 1 so that only non-lake cells are 0
    arcpy.AddField_management(zone_path, zone_value_field, 'LONG')
    arcpy.AddField_management(zone_path, int_area_field, 'DOUBLE')
    arcpy.AddField_management(zone_path, int_pct_field, 'DOUBLE')
    if zone_value == arcpy.Describe(zone_path).OIDFieldName:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}! + 1'.format(zone_value), 'PYTHON')

    # If zone value is an INT, save it into a field first
    elif type(zone_value) is int:
        # zone_value = int(zone_value)
        arcpy.CalculateField_management(
            zone_path, zone_value_field, zone_value, 'PYTHON')
    # Use zone_value field directly
    else:
        # zone_value_field = zone_value
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}!'.format(zone_value), 'PYTHON')

    # Calculate area of HRU cell if necessary
    # if not arcpy.ListFields(zone_path, area_field):
    #    arcpy.AddField_management(zone_path, area_field, 'DOUBLE')
    #    cell_area_func(zone_path, area_field)

    # Intersect the zone layer with the HRU
    # zone_int_path = os.path.join('in_memory', 'hru_lakes')
    zone_int_path = zone_path.replace('.shp', '_intersect.shp')
    arcpy.Intersect_analysis(
        (hru_param_path, zone_path), zone_int_path)

    # Calculate using cell_area_func to force units to match
    cell_area_func(zone_int_path, int_area_field)

    n = int(arcpy.GetCount_management(zone_int_path).getOutput(0))
    block_size = 200000
    for i, x in enumerate(range(0, n, block_size)):
        logging.debug('    FIDS: {0}-{1}'.format(x, x + block_size))
        subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
            hru_param.fid_field, x, x + block_size)
            # arcpy.Describe(zone_int_path).OIDFieldName, x, x + block_size)

        # Read in FID of selected cells
        hru_cell_dict = dict()
        fields = [
            hru_param.fid_field, hru_area_field,
            int_area_field, zone_value_field]
        with arcpy.da.SearchCursor(zone_int_path, fields, subset_str) as s_cursor:
            for row in s_cursor:
                if (100 * float(row[2]) / float(row[1])) >= area_pct:
                    hru_cell_dict[int(row[0])] = [float(row[3]), float(row[2])]

        # Set value of selected HRU cells
        fields = [hru_param.fid_field, zone_field]
        if zone_area_field:
            fields.append(zone_area_field)
        with arcpy.da.UpdateCursor(hru_param_path, fields, subset_str) as u_cursor:
            for row in u_cursor:
                # Remove items to speed up subsequent searches
                try:
                    if len(fields) == 3:
                        row[1], row[2] = hru_cell_dict.pop(int(row[0]))
                    elif len(fields) == 2:
                        row[1] = hru_cell_dict.pop(int(row[0]))
                    u_cursor.updateRow(row)
                except KeyError:
                    pass
        del hru_cell_dict


def zone_by_centroid_func(zone_path, zone_field, zone_value,
                          hru_param_path, hru_point_path, hru_param):
    """Flag cells that are inside a feature based on the centroid location

    Set values that are in zone, but don't reset values that are out of zone

    Args:
        zone_path (str):
        zone_field (str):
        zone_value (int):
        hru_param_path (str):
        hru_point_path (str):
        hru_param: class:`support_functions.HRUParameters`

    Returns:
        None
    """
    logging.debug('\nzone_by_centroid_func')
    logging.debug('  {}'.format(zone_path))
    # Need to set zone value into a field before intersect
    # If zone_value is FID, add 1 so that only zone cells are 0
    zone_value_field = 'ZONE_VALUE'
    arcpy.AddField_management(zone_path, zone_value_field, 'LONG')
    if zone_value == arcpy.Describe(zone_path).OIDFieldName:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}! + 1'.format(zone_value), 'PYTHON')
    # Save zone value into a field first
    elif type(zone_value) is int:
        arcpy.CalculateField_management(
            zone_path, zone_value_field, zone_value, 'PYTHON')
        # DEADBEEF
        # Use zone_value field directly
        # zone_value = int(zone_value)
    else:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}!'.format(zone_value), 'PYTHON')
        # DEADBEEF
        # Use zone_value field directly
        # zone_value_field = zone_value

    # Intersect the zone layer with the HRU
    # zone_int_path = os.path.join('in_memory', 'hru_ppt_zones')
    zone_int_path = zone_path.replace('.shp', '_intersect.shp')
    arcpy.Intersect_analysis(
        (hru_point_path, zone_path), zone_int_path)

    # DEADBEEF - Why do I make a layer and select all features?
    zone_int_layer = 'zone_int_layer'
    arcpy.MakeFeatureLayer_management(zone_int_path, zone_int_layer)
    arcpy.SelectLayerByAttribute_management(zone_int_layer, 'CLEAR_SELECTION')
    arcpy.SelectLayerByAttribute_management(zone_int_layer, 'SWITCH_SELECTION')

    n = int(arcpy.GetCount_management(zone_int_layer).getOutput(0))
    block_size = 200000
    for i, x in enumerate(range(0, n, block_size)):
        logging.debug('    FIDS: {0}-{1}'.format(x, x + block_size))
        subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
            hru_param.fid_field, x, x + block_size)
            # arcpy.Describe(zone_int_layer).OIDFieldName, x, x + block_size)

        # Read in FID of selected cells
        hru_point_dict = dict()
        fields = (hru_param.fid_field, zone_value_field)
        with arcpy.da.SearchCursor(zone_int_layer, fields, subset_str) as s_cursor:
            for row in s_cursor:
                hru_point_dict[int(row[0])] = row[1]

        # Set value of selected HRU cells
        fields = (hru_param.fid_field, zone_field)
        with arcpy.da.UpdateCursor(hru_param_path, fields, subset_str) as u_cursor:
            for row in u_cursor:
                # Remove items to speed up subsequent searches
                try:
                    row[1] = hru_point_dict.pop(int(row[0]))
                    u_cursor.updateRow(row)
                except KeyError:
                    pass
        del hru_point_dict

    # Cleanup
    arcpy.Delete_management(zone_int_layer)


def jensen_haise_func(hru_param_path, jh_coef_field, hru_elev_field,
                      jh_tmin_field, jh_tmax_field):
    """"""
    jh_cb = (
        'def ea(temp_c):\n' +
        '    return 6.1078 * math.exp((17.269 * temp_c) / (temp_c + 237.3))\n' +
        'def jensen_haise(elev, t_low, t_high):\n' +
        '    return 27.5 - 0.25 * (ea(t_high) - ea(t_low)) - (elev / 1000)\n')
    arcpy.CalculateField_management(
        hru_param_path, jh_coef_field,
        'jensen_haise(!{0}!, !{1}!, !{2}!)'.format(
            hru_elev_field, jh_tmin_field, jh_tmax_field),
        'PYTHON', jh_cb)


def remap_check(remap_path):
    """"""
    # Check that the file exists
    if not os.path.isfile(remap_path):
        logging.error(
            '\nERROR: ASCII remap file ({0}) does not exist\n'.format(
                os.path.basename(remap_path)))
        sys.exit()

    # Read in the remap
    with open(remap_path, 'r') as remap_f:
        remap_lines = remap_f.readlines()
        line_count = len(remap_lines)

    # Problems reading/applying ASCII remap files can be caused by:
    #   Blank lines
    #   Empty lines at the end (from a final newline character)
    #   ArcGIS 10.2 - Old style in line comments (/*)
    #   ArcGIS 10.2 - Comments can't be in first line (?)
    #   ArcGIS 10.2 - Comments can't be longer than 80 characters (?)
    # If either of these are present in the file, resave the filtered lines

    # Check for old style comments (/*) in ASCII remap files
    # This could be changed to save the comments at the end of the file
    if arcpy.GetInstallInfo()['Version'].startswith('10.2'):
        if any([l for l in remap_lines if '/*' in l]):
            logging.error(
                    ('\nERROR: ASCII remap file ({0}) has pre-ArcGIS 10.2 ' +
                     'comments (\*)\n  Try running the "convert_remap_arc10p2.py"' +
                     'script\n').format(os.path.basename(remap_path)))
            sys.exit()

    # First check for final newline character
    save_flag = False
    if remap_lines and remap_lines[-1] and remap_lines[-1].endswith('\n'):
        logging.debug('  Final newline character')
        save_flag = True

    # Then remove empty lines and strip white space and newline characters
    # If lines were removed, resave the filtered remap file
    remap_lines = [l.strip() for l in remap_lines]
    remap_lines = [l for l in remap_lines if l]
    if len(remap_lines) != line_count:
        logging.debug('  Whitespace or empty lines')
        save_flag = True

    # Trim comments longer than 80 characters
    if arcpy.GetInstallInfo()['Version'].startswith('10.2'):
        if any([len(l) > 80 for l in remap_lines if "#" in l]):
            remap_lines = [l[:79] if "#" in l else l for l in remap_lines]
            logging.debug('  Lines longer than 80 characters')
            save_flag = True

    # If lines were removed, resave the filtered remap file
    if save_flag:
        logging.warning(
            '  The ASCII remap file ({0}) will be overwritten'.format(
                os.path.basename(remap_path)))
        with open(remap_path, 'w') as remap_f:
            for i, line in enumerate(remap_lines):
                # Don't write newline character on last line
                # This causes an error in ArcGIS 10.2.2
                if (i+1) < len(remap_lines):
                    remap_f.write(line + '\n')
                else:
                    remap_f.write(line)
    return True


#  Remap aspect
# logging.info('\nRemapping Aspect to HRU_ASPECT')
# arcpy.CalculateField_management(
#    polygon_path, hru_aspect_field,
#    'Reclass(!{0}!)'.format(dem_aspect_field),
#    'PYTHON', remap_code_block(aspect_remap_path))
# # arcpy.DeleteField_management(polygon_path, dem_aspect_field)


def remap_code_block(remap_path):
    """"""
    with open(remap_path) as remap_f:
        lines = remap_f.readlines()
    remap_cb = ''
    for l in lines:
        # Skip comment lines
        if '#' in l:
            continue
        # Remove remap description
        l = l.strip().split('/*')[0]
        # Split line on spaces and semi-colon
        l_split = [item.strip() for item in re.split('[ :]+', l)]
        # Remap as a range if a min, max and value are all present
        if len(l_split) == 3:
            range_remap_flag = True
        # Otherwise remap directly
        elif len(l_split) == 2:
            range_remap_flag = False
        # Skip lines that don't match format
        else:
            continue
        # Write remap code block
        if not range_remap_flag:
            if not remap_cb:
                remap_cb = ('    if value == {0}: ' +
                            'return {1}\n'.format(*l_split))
            else:
                remap_cb += ('    elif value == {0}: ' +
                             'return {1}\n'.format(*l_split))
        else:
            if not remap_cb:
                remap_cb = ('    if (value >= {0} and value <= {1}): ' +
                            'return {2}\n').format(*l_split)
            else:
                remap_cb += ('    elif (value > {0} and value <= {1}): ' +
                             'return {2}\n').format(*l_split)
    remap_cb = 'def Reclass(value):\n' + remap_cb
    return remap_cb


# def reclass_ascii_float_func(raster_path, remap_path):
#    # Read remap file into memory
#    with open(remap_path) as remap_f: lines = remap_f.readlines()
#    remap_f.close()
#    first_line = True
#    raster_obj = Raster(raster_path)
#    for l in lines:
#        # Skip comment lines
#        if '#' in l: continue
#        # Remove remap description
#        l = l.split('/*')[0]
#        # Split line on spaces and semi-colon
#        l_split = map(float, [item for item in re.split('[ :]+', l) if item])
#        # Remap as a range if a min, max and value are all present
#        if len(l_split) == 3: range_remap_flag = True
#        # Otherwise remap directly
#        elif len(l_split) == 2: range_remap_flag = False
#        # Skip lines that don't match format
#        else: continue
#        # Write remap code block
#        if not range_remap_flag:
#            raster_obj = Con(raster_obj == l_split[0], l_split[1], raster_obj)
#        elif first_line:
#            raster_obj = Con(
#                ((raster_obj >= l_split[0]) & (raster_obj <= l_split[1])),
#                l_split[2], raster_obj)
#        else:
#            raster_obj = Con(
#                ((raster_obj > l_split[0]) & (raster_obj <= l_split[1])),
#                l_split[2], raster_obj)
#        first_line = False
#    return raster_obj


def is_number(s):
    """"""
    try:
        float(s)
        return True
    except ValueError:
        return False


def raster_path_to_array(input_path, mask_extent=None, return_nodata=False):
    """"""
    return raster_obj_to_array(Raster(input_path), mask_extent, return_nodata)

def raster_obj_to_array(input_obj, mask_extent=None, return_nodata=False):
    """"""
    input_nodata = input_obj.noDataValue
    input_cs = input_obj.meanCellHeight
    input_rows, input_cols = input_obj.height, input_obj.width
    input_extent = input_obj.extent
    if mask_extent:
        int_extent = get_extent_intersection([input_extent, mask_extent])
        int_pnt = arcpy.Point()
        int_pnt.X = int_extent.XMin
        int_pnt.Y = int_extent.YMin
        int_rows, int_cols = extent_shape(int_extent, input_cs)
        output_array = arcpy.RasterToNumPyArray(
            input_obj, int_pnt, int_cols, int_rows)
    else:
        output_array = arcpy.RasterToNumPyArray(input_obj)
    # Integer type raster can't have NaN values, will only set floats to NaN
    if (output_array.dtype == np.float32 or
        output_array.dtype == np.float64):
        output_array[output_array == input_nodata] = np.NaN
        output_nodata = np.NaN
    else:
        output_nodata = int(input_nodata)
    if return_nodata:
        return output_array, output_nodata
    else:
        return output_array


def array_to_raster(input_array, output_path, pnt, cs, mask_array=None):
    """"""
    output_array = np.copy(input_array)
    # Float arrays have to have nodata set to some value (-9999)
    if (output_array.dtype == np.float32 or
        output_array.dtype == np.float64):
        output_nodata = -9999
        output_array[np.isnan(output_array)] = output_nodata
    # Boolean arrays need to be converted unsigned ints
    elif output_array.dtype == np.bool:
        output_array = output_array.astype(np.uint8)
        output_nodata = 255
    # If a mask array is give, assume all 0 values are nodata
    if np.any(mask_array): output_array[mask_array == 0] = output_nodata
    output_obj = arcpy.NumPyArrayToRaster(
        output_array, pnt, cs, cs, output_nodata)
    output_obj.save(output_path)
    del output_obj
    arcpy.DefineProjection_management(
        output_path, env.outputCoordinateSystem)
    arcpy.CalculateStatistics_management(output_path)


def flood_fill(test_array, four_way_flag=True, edge_flt=None):
    """Flood fill algorithm"""
    input_array = np.copy(test_array)
    input_rows, input_cols = input_array.shape
    h_max = np.nanmax(input_array * 2.0)
    logging.debug("  Hmax: %s" % (h_max/2.0))

    # Since ArcGIS doesn't ship with SciPy (only numpy), don't use ndimage module
    if four_way_flag:
        el = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]).astype(np.bool)
    else:
        el = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool)
    # if four_way_flag:
    #    el = ndimage.generate_binary_structure(2,1).astype(np.int)
    # else:
    #    el = ndimage.generate_binary_structure(2,2).astype(np.int)

    # Build data/inside/edge masks
    data_mask = ~np.isnan(input_array)

    # Since ArcGIS doesn't ship with SciPy (only numpy), don't use ndimage module
    inside_mask = np_binary_erosion(data_mask, structure=el)
    # inside_mask = ndimage.binary_erosion(data_mask, structure=el)

    edge_mask = (data_mask & ~inside_mask)
    # Initialize output array as max value test_array except edges
    output_array = np.copy(input_array)
    output_array[inside_mask] = h_max
    # Set edge pixels less than edge_flt to edge_flt
    if edge_flt:
        output_array[edge_mask & (output_array<=edge_flt)]=edge_flt

    # Build priority queue and place edge pixels into queue
    put = heapq.heappush
    get = heapq.heappop
    fill_heap = [
        (output_array[t_row,t_col], int(t_row), int(t_col), 1)
        for t_row, t_col in np.transpose(np.where(edge_mask))]
    heapq.heapify(fill_heap)
    # logging.info("    Queue Size: %s" % len(fill_heap))

    # Cleanup
    del data_mask, edge_mask, el
    # logging.info("    Prep Time: %s" % (clock()-start_total))

    while True:
        try:
            h_crt, t_row, t_col, edge_flag = get(fill_heap)
        except IndexError:
            break
        for n_row, n_col in [
            ((t_row-1), t_col), ((t_row+1), t_col),
            (t_row, (t_col-1)), (t_row, (t_col+1))]:
            # Skip cell if outside array edges
            if edge_flag:
                try:
                    if not inside_mask[n_row, n_col]:
                        continue
                except IndexError:
                    continue
            if output_array[n_row, n_col]==h_max:
                output_array[n_row, n_col] = max(
                    h_crt, input_array[n_row, n_col])
                put(fill_heap,
                    (output_array[n_row, n_col], n_row, n_col, 0))
    return output_array


def np_binary_erosion(input_array, structure=np.ones((3,3)).astype(np.bool)):
    """NumPy binary erosion function

    No error checking on input array (type)
    No error checking on structure element (# of dimensions, shape, type, etc.)

    Args:
        input_array: Binary NumPy array to be eroded. Non-zero (True) elements
            form the subset to be eroded
        structure: Structuring element used for the erosion. Non-zero elements
            are considered True. If no structuring element is provided, an
            element is generated with a square connectivity equal to one.

    Returns:
        binary_erosion: Erosion of the input by the stucturing element
    """
    rows, cols = input_array.shape

    # Pad output array (binary_erosion) with extra cells around the edge
    # so that structuring element will fit without wrapping.
    # A 3x3 structure, will need 1 additional cell around the edge
    # A 5x5 structure, will need 2 additional cells around the edge
    output_shape = tuple(
        ss + dd - 1 for ss, dd in zip(input_array.shape, structure.shape))
    input_pad_array = np.zeros(output_shape).astype(np.bool)
    input_pad_array[1:rows+1, 1:cols+1] = input_array
    binary_erosion = np.zeros(output_shape).astype(np.bool)

    # Cast structure element to boolean
    struc_mask = structure.astype(np.bool)

    # Iterate over each cell
    for row in xrange(rows):
        for col in xrange(cols):
            # The value of the output pixel is the minimum value of all the
            #   pixels in the input pixel's neighborhood.
            binary_erosion[row+1, col+1] = np.min(
                input_pad_array[row:row+3, col:col+3][struc_mask])
    return binary_erosion[1:rows+1, 1:cols+1]
