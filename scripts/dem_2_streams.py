#--------------------------------
# Name:         dem_2_streams.py
# Purpose:      GSFLOW Flow Parameters
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2015-06-05
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import heapq
import logging
import os
import re
import sys
from time import clock

import arcpy
from arcpy import env
from arcpy.sa import *

import numpy as np
##from scipy import ndimage

from support_functions import *

################################################################################

def gsflow_flow_parameters(config_path, overwrite_flag=False, debug_flag=False):
    """Calculate GSFLOW Flow Parameters

    Args:
        config_file: Project config file path
        ovewrite_flag: boolean, overwrite existing files
        debug_flag: boolean, enable debug level logging
    Returns:
        None
    """

    try:
        logging.info('\nGSFLOW Flow Parameters')

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
        logging.debug('\nReading Input File')

        ## Log DEBUG to file
        log_file_name = 'gsflow_dem_2_stream_log.txt'
        log_console = logging.FileHandler(
            filename=os.path.join(hru.log_ws, log_file_name), mode='w')
        log_console.setLevel(logging.DEBUG)
        log_console.setFormatter(logging.Formatter('%(message)s'))
        logging.getLogger('').addHandler(log_console)
        logging.info('\nGSFLOW DEM To Streams')

        ## Check whether lake parameters should be calculated
        try:
            set_lake_flag = inputs_cfg.getboolean('INPUTS', 'set_lake_flag')
        except:
            logging.debug('  set_lake_flag = False')
            set_lake_flag = False

        ## Check whether ocean parameters should be calculated
        try:
            set_ocean_flag = inputs_cfg.getboolean('INPUTS', 'set_ocean_flag')
        except:
            logging.debug('  set_ocean_flag = False')
            set_ocean_flag = False

        ## Subbasin points
        subbasin_input_path = inputs_cfg.get('INPUTS', 'subbasin_points_path')
        subbasin_zone_field = inputs_cfg.get('INPUTS', 'subbasin_zone_field')

        ## Flow parameters
        reset_dem_adj_flag = inputs_cfg.getboolean('INPUTS', 'reset_dem_adj_flag')
        calc_flow_dir_points_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_flow_dir_points_flag')
        calc_sinks_8_way_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_sinks_8_way_flag')
        calc_sinks_4_way_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_sinks_4_way_flag')
        flow_acc_threshold = inputs_cfg.getint('INPUTS', 'flow_acc_threshold')
        flow_length_threshold = inputs_cfg.getint('INPUTS', 'flow_length_threshold')
        try: 
    	    lake_seg_offset = inputs_cfg.getint('INPUTS', 'lake_seg_offset')
        except: 
            lake_seg_offset = 0
        if lake_seg_offset < 0:
            logging.error(
                '\nERROR: lake_seg_offset must be an integer greater than 0')
            raise SystemExit()

        ## Check input paths
        dem_temp_ws = os.path.join(hru.param_ws, 'dem_rasters')
        dem_path = os.path.join(dem_temp_ws, 'dem.img')
        if not arcpy.Exists(dem_path):
            logging.error(
                ('\nERROR: Projected/clipped DEM ({0}) does not exist'+
                '\nERROR: Try rerunning gsflow_dem_parameters.py').format(dem_path))
            raise SystemExit()
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(hru.polygon_path))
            raise SystemExit()

        ## Check subbasin_points
        logging.info('\nChecking input subbasin shapefile')
        if not os.path.isfile(subbasin_input_path): 
            logging.error(
                ('\nERROR: Subbasin points shapefiles does not exist'+
                 '\nERROR:   {0}').format(subbasin_input_path))
            raise SystemExit()
        ## subbasin_zone_path must be a point shapefile
        elif arcpy.Describe(subbasin_input_path).datasetType <> 'FeatureClass':
            logging.error(
                '\nERROR: subbasin_input_path must be a point shapefile')
            raise SystemExit()
        ## Check subbasin_zone_fields
        if subbasin_zone_field.upper() in ['', 'FID', 'NONE']:
            subbasin_fid_field = arcpy.Describe(subbasin_input_path).OIDFieldName
            logging.warning(
                '  NOTE: Using {0}+1 to set {1}\n'.format(
                    subbasin_fid_field, hru.subbasin_field))
            subbasin_zone_field = 'ZONE_VALUE'
            if not arcpy.ListFields(subbasin_input_path, subbasin_zone_field):
                arcpy.AddField_management(subbasin_input_path, subbasin_zone_field, 'LONG')
            arcpy.CalculateField_management(
                subbasin_input_path, subbasin_zone_field,
                '!{0}! + 1'.format(subbasin_fid_field), 'PYTHON')
        elif not arcpy.ListFields(subbasin_input_path, subbasin_zone_field):
            logging.error(
                '\nERROR: subbasin_zone_field {0} does not exist\n'.format(
                    subbasin_zone_field))
            raise SystemExit()
        ## Need to check that subbasin_zone_field is an int type
        elif not [f.type for f in arcpy.Describe(subbasin_input_path).fields
                  if (f.name == subbasin_zone_field and
                      f.type in ['SmallInteger', 'Integer'])]:
            logging.error(
                '\nERROR: subbasin_zone_field {0} must be an integer type\n'.format(
                    subbasin_zone_field))
            raise SystemExit()
        ## Need to check that subbasin_zone_field is all positive values
        if min([row[0] for row in arcpy.da.SearchCursor(
            subbasin_input_path, [subbasin_zone_field])]) <= 0:
            logging.error(
                '\nERROR: subbasin_zone_field values must be positive\n'.format(
                    subbasin_zone_field))
            raise SystemExit()

        ## Build output folder if necessary
        flow_temp_ws = os.path.join(hru.param_ws, 'flow_rasters')
        if not os.path.isdir(flow_temp_ws):
            os.mkdir(flow_temp_ws)
        ## Output paths
        hru_type_in_path = os.path.join(flow_temp_ws, 'hru_type_in.img')
        hru_type_path = os.path.join(flow_temp_ws, 'hru_type.img')
        dem_adj_path = os.path.join(flow_temp_ws, 'dem_adj.img')
        lake_id_path = os.path.join(flow_temp_ws, 'lake_id.img')
        dem_sink8_path = os.path.join(flow_temp_ws, 'dem_sink_8_way.img')
        dem_sink4_path = os.path.join(flow_temp_ws, 'dem_sink_4_way.img')
        dem_fill_path = os.path.join(flow_temp_ws, 'dem_fill.img')
        flow_dir_path = os.path.join(flow_temp_ws, 'flow_dir.img')
        flow_dir_points = os.path.join(flow_temp_ws, 'flow_dir_points.shp')
        flow_acc_full_path = os.path.join(flow_temp_ws, 'flow_acc_full.img')
        flow_acc_sub_path = os.path.join(flow_temp_ws, 'flow_acc_sub.img')
        flow_mask_path = os.path.join(flow_temp_ws, 'flow_mask.img')
        stream_link_path = os.path.join(flow_temp_ws, 'stream_link.img')
        stream_link_a_path = os.path.join(flow_temp_ws, 'stream_link_a.img')
        stream_link_b_path = os.path.join(flow_temp_ws, 'stream_link_b.img')
        stream_order_path = os.path.join(flow_temp_ws, 'stream_order.img')
        stream_length_path = os.path.join(flow_temp_ws, 'stream_length.img')
        watersheds_path = os.path.join(flow_temp_ws, 'watersheds.img')
        subbasin_points_path = os.path.join(flow_temp_ws, 'subbasin_points.shp')
        subbasin_path = os.path.join(flow_temp_ws, 'subbasin.img')
        basin_path = os.path.join(flow_temp_ws, 'basin.img')
        streams_path = os.path.join(flow_temp_ws, 'streams.shp')

        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        ##env.pyramid = 'PYRAMIDS -1'
        env.pyramid = 'PYRAMIDS 0'
        env.workspace = flow_temp_ws
        env.scratchWorkspace = hru.scratch_ws

        ## Check DEM field
        logging.info('\nAdding DEM fields if necessary')
        add_field_func(hru.polygon_path, hru.iseg_field, 'LONG')
        add_field_func(hru.polygon_path, hru.irunbound_field, 'LONG')
        add_field_func(hru.polygon_path, hru.flow_dir_field, 'LONG')
        add_field_func(hru.polygon_path, hru.dem_sink8_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_sink4_field, 'DOUBLE')

        ## Set environment parameters
        env.extent = hru.extent
        env.cellsize = hru.cs
        env.outputCoordinateSystem = hru.sr
        
        ## Check that dem_adj_copy_field exists
        ##if len(arcpy.ListFields(hru.polygon_path, dem_adj_copy_field)) == 0:
        ##    logging.error('\nERROR: dem_adj_copy_field {0} does not exist\n'.format(
        ##        dem_adj_copy_field))
        ##    raise SystemExit()
        ## Reset DEM_ADJ
        ##if reset_dem_adj_flag:
        ##    logging.info('\nResetting {0} to {1}'.format(
        ##        dem_adj_field, dem_adj_copy_field))
        ##    arcpy.CalculateField_management(
        ##        hru.polygon_path, dem_adj_field,
        ##        '!{0}!'.format(dem_adj_copy_field), 'PYTHON')


        ## Check lake cell elevations
        if set_lake_flag:
            logging.info('\nChecking lake cell {0}'.format(hru.dem_adj_field))
            lake_elev_dict = defaultdict(list)
            fields = [
                hru.type_in_field, hru.lake_id_field,
                hru.dem_adj_field, hru.id_field]
            for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
                if int(row[0]) <> 2:
                    continue
                lake_elev_dict[int(row[1])].append(float(row[2]))
            del fields
            logging.info('  {0:>7} {1:>12} {2:>12} {3:>12} {4:>12}'.format(
                'Lake ID', 'Minimum', 'Mean', 'Maximum', 'Std. Dev.'))
            for lake_id, lake_elev_list in lake_elev_dict.items():
                lake_elev_array = np.array(lake_elev_list)
                logging.info('  {0:7} {1:12f} {2:12f} {3:12f} {4:12f}'.format(
                    lake_id, np.min(lake_elev_array), np.mean(lake_elev_array),
                    np.max(lake_elev_array), np.std(lake_elev_array)))
                if np.std(lake_elev_array) > 1:
                    logging.warning(
                        '  Please check the lake cell elevations\n'+
                        '  They may need to be manually adjusted'.format(lake_id))
                    raw_input('  Press ENTER to continue')
                del lake_elev_array


        ## Convert DEM_ADJ to raster
        logging.info('\nExporting HRU polygon parameters to raster')
        logging.debug('  HRU_TYPE')
        ## Read in original HRU_TYPE for defining subbasins
        arcpy.PolygonToRaster_conversion(
            hru.polygon_path, hru.type_in_field, hru_type_in_path,
            "CELL_CENTER", "", hru.cs)
        hru_type_in_obj = Raster(hru_type_in_path)
        logging.debug('  DEM_ADJ')
        arcpy.PolygonToRaster_conversion(
            hru.polygon_path, hru.dem_adj_field, dem_adj_path,
            "CELL_CENTER", "", hru.cs)
        dem_adj_obj = Raster(dem_adj_path)
        if set_lake_flag:
            logging.debug('  LAKE_ID')
            arcpy.PolygonToRaster_conversion(
                hru.polygon_path, hru.lake_id_field, lake_id_path,
                "CELL_CENTER", "", hru.cs)
            lake_id_obj = Raster(lake_id_path)

        ## Flow Direction
        logging.info('\nCalculating flow direction')
        flow_dir_obj = FlowDirection(dem_adj_obj, True)
        ##flow_dir_obj = FlowDirection(dem_adj_obj, False)
        flow_dir_obj.save(flow_dir_path)
        

        ## Fill DEM_ADJ raster
        logging.info('Filling DEM_ADJ (8-way)')
        dem_fill_obj = Fill(dem_adj_obj)
        dem_fill_obj.save(dem_fill_path)

        ## Sinks (8-way)
        if calc_sinks_8_way_flag:
            logging.info('Calculating sinks (8-way)')
            dem_sink8_obj = Con(
                ~IsNull(Sink(flow_dir_obj)), (dem_fill_obj - dem_adj_obj))
            dem_sink8_obj.save(dem_sink8_path)
            dem_sink8_array = raster_to_array(dem_sink8_path)[0]
            if np.all(np.isnan(dem_sink8_array)):
                logging.info('  No sinks (8-way)')
                fill_flag = False
            else:
                fill_flag = True
            del dem_sink8_array, dem_sink8_obj, dem_adj_obj
        
        ## Sinks (4-way)
        if calc_sinks_4_way_flag:
            logging.info('Calculating sinks (4-way)')
            pnt = arcpy.Point()
            pnt.X = hru.extent.XMin
            pnt.Y = hru.extent.YMin
            dem_adj_array = raster_to_array(dem_adj_path)[0].astype(np.float)
            dem_sink4_array = flood_fill(dem_adj_array, True)
            dem_sink4_array -= dem_adj_array
            dem_sink4_array[dem_sink4_array == 0] = np.nan
            array_to_raster(dem_sink4_array, dem_sink4_path, pnt, hru.cs)
            if np.all(np.isnan(dem_sink4_array)):
                logging.info('  No sinks (4-way)')
                fill_flag = False
            else:
                fill_flag = True
            del dem_adj_array, dem_sink4_array, pnt

        ## Recaculate Flow Direction
        if fill_flag:
            logging.info('Re-calculating flow direction')
            flow_dir_obj = FlowDirection(dem_fill_obj, True)
            flow_dir_obj.save(flow_dir_path)

        ## Save flow direction as points
        if calc_flow_dir_points_flag:
            logging.info('Flow direction points')
            ## ArcGIS fails for raster_to_x conversions on a network path
            ## You have to go through an in_memory file first
            flow_dir_temp = os.path.join('in_memory', 'flow_dir')
            arcpy.RasterToPoint_conversion(flow_dir_obj, flow_dir_temp)
            arcpy.CopyFeatures_management(flow_dir_temp, flow_dir_points)
            arcpy.Delete_management(flow_dir_temp)
            del flow_dir_temp
            ## Reclassify flow directions to angles, assuming 1 is 0
            remap_cb = (
                'def Reclass(value):\n'+
                '    if value == 1: return 0\n'+
                '    elif value == 2: return 45\n'+
                '    elif value == 4: return 90\n'+
                '    elif value == 8: return 135\n'+
                '    elif value == 16: return 180\n'+
                '    elif value == 32: return 225\n'+
                '    elif value == 64: return 270\n'+
                '    elif value == 128: return 315\n')
            arcpy.CalculateField_management(
                flow_dir_points, 'grid_code',
                'Reclass(!{0}!)'.format('grid_code'), 'PYTHON', remap_cb)

        ## Write flow direction to hru_polygon
        logging.debug('  Extracting flow direction at points')
        vt_list = [[flow_dir_path, hru.flow_dir_field]]
        mem_point_path = os.path.join('in_memory', 'hru_point')
        arcpy.CopyFeatures_management(hru.point_path, mem_point_path)
        ExtractMultiValuesToPoints(mem_point_path, vt_list, 'NONE')
        logging.debug('  Reading flow direction values at point')
        data_dict = defaultdict(dict)
        fields = [hru.flow_dir_field, hru.fid_field]
        with arcpy.da.SearchCursor(mem_point_path, fields) as s_cursor:
            for row in s_cursor:
                ## Set nodata cells to 0
                if row[1] is not None:
                    data_dict[int(row[1])][hru.flow_dir_field] = int(row[0])
                del row
        logging.debug('  Writing flow direction values to polygon')
        fields = [hru.flow_dir_field, hru.fid_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                row_dict = data_dict.get(int(row[-1]), None)
                for i, field in enumerate(fields[:-1]):
                    if row_dict:
                        row[i] = row_dict[field]
                    else:
                        row[i] = 0   
                u_cursor.updateRow(row)
                del row_dict, row



        ## Subbasins
        logging.info('\nChecking input subbasin points')
        ## Check that subbasin values increment from 1 to nsub
        logging.debug('  Checking subbasin ID')
        subbasin_id_list = sorted(list(set(
            [row[0] for row in arcpy.da.SearchCursor(
                subbasin_input_path, [subbasin_zone_field])])))
        if subbasin_id_list <> range(1,len(subbasin_id_list)+1):
            logging.error(
                ('\nERROR: SUB_BASINs must be sequential starting from 1'+
                 '\nERROR:   {0}').format(subbasin_id_list))
            raise SystemExit()
        subbasin_input_count = len(subbasin_id_list)
        logging.debug('  {} subbasins'.format(subbasin_input_count))
        ## Get spatial reference of subbasin_points
        subbasin_points_desc = arcpy.Describe(subbasin_input_path)
        subbasin_points_sr = subbasin_points_desc.spatialReference
        logging.debug('  Subbasin points spat. ref.:  {0}'.format(
            subbasin_points_sr.name))
        logging.debug('  Subbasin points GCS:         {0}'.format(
            subbasin_points_sr.GCS.name))
        if arcpy.Exists(subbasin_points_path):
            arcpy.Delete_management(subbasin_points_path)
        ## Project points if necessary
        if hru.sr.name <> subbasin_points_sr.name:
            ## Set preferred transforms
            transform_str = transform_func(hru.sr, subbasin_points_sr)
            logging.debug('    Transform: {0}'.format(transform_str))
            ## Project subbasin_points to match HRU_Polygon
            logging.debug('  Projecting subbasin_points')
            logging.debug('    {0}'.format(subbasin_input_path))
            logging.debug('    {0}'.format(subbasin_points_path))
            arcpy.ClearEnvironment("outputCoordinateSystem")
            arcpy.ClearEnvironment("extent")
            ##env.scratchWorkspace = scratch_ws
            ##arcpy.ClearEnvironment("cellsize")
            arcpy.Project_management(
                subbasin_input_path, subbasin_points_path, hru.sr,
                transform_str, subbasin_points_sr)
            env.extent = hru.extent
            env.outputCoordinateSystem = hru.sr
            ##env.scratchWorkspace = 'in_memory'
            ##env.cellsize = hru.cs
        else:
            arcpy.Copy_management(subbasin_input_path, subbasin_points_path)
            
        ## Select the HRU cells that intersect the subbasin point cells
        logging.debug('  Reading input subbasin points')
        hru_polygon_lyr = "hru_polygon_lyt"
        arcpy.MakeFeatureLayer_management(hru.polygon_path, hru_polygon_lyr)
        arcpy.SelectLayerByLocation_management (
            hru_polygon_lyr, "intersect", subbasin_points_path)
        input_xy_dict = dict()
        fields = [hru.col_field, hru.row_field, hru.x_field, hru.y_field]
        for row in arcpy.da.SearchCursor(hru_polygon_lyr, fields):
            cell = (int(row[0]), int(row[1]))
            input_xy_dict[cell] = (int(row[2]), int(row[3]))
        arcpy.Delete_management(hru_polygon_lyr)
        del hru_polygon_lyr
        ##for k,v in input_xy_dict.items():
        ##    logging.debug('    {0} {1}'.format(k,v))
        
        ## First calculate downstream cell for all cells
        logging.info('\nBuilding all subbasin points')
        out_cell_dict = dict()
        hru_type_in_dict = dict()
        cell_xy_dict = dict()
        fields = [
            hru.type_in_field, hru.flow_dir_field, hru.id_field,
            hru.col_field, hru.row_field, hru.x_field, hru.y_field]
        for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
            cell = (int(row[3]), int(row[4]))
            out_cell_dict[cell] = next_row_col(int(row[1]), cell)
            hru_type_in_dict[cell] = int(row[0])
            cell_xy_dict[cell] = (int(row[5]), int(row[6]))

        ## Identify all active/lake cells that flow to inactive water cells
        exit4_xy_list = sorted([
            cell_xy for cell, cell_xy in cell_xy_dict.items()
            if (cell not in input_xy_dict.keys() and
                cell in hru_type_in_dict.keys() and
                hru_type_in_dict[cell] in [1,2] and
                hru_type_in_dict[out_cell_dict[cell]] == 4)])
        ## Identify all active/lake cells that flow to inactive land cells
        exit0_xy_list = sorted([
            cell_xy for cell, cell_xy in cell_xy_dict.items()
            if (cell not in input_xy_dict.keys() and
                cell in hru_type_in_dict.keys() and
                hru_type_in_dict[cell] in [1,2] and
                hru_type_in_dict[out_cell_dict[cell]] == 0)])
        ## DEADBEEF - Separate subbasin for active/lake cells that flow out of the grid?
        ##exit_xy_list = sorted([
        ##    cell_xy for cell, cell_xy in cell_xy_dict.items()
        ##    if (cell not in input_xy_dict.keys() and
        ##        cell in hru_type_in_dict.keys() and hru_type_in_dict[cell] in [1,2] and
        ##        out_cell_dict[cell] not in hru_type_in_dict.keys())])
        ## Create subbasin points for cells
        fields = ["SHAPE@XY", subbasin_zone_field]
        with arcpy.da.InsertCursor(subbasin_points_path, fields) as insert_c:
            for out_cell_xy in exit4_xy_list:
                insert_c.insertRow([out_cell_xy, subbasin_input_count+1])
            for out_cell_xy in exit0_xy_list:
                insert_c.insertRow([out_cell_xy, subbasin_input_count+2])
                ##out_cell_xy, out_cell_i+subbasin_input_count+2])
        del exit4_xy_list, exit0_xy_list, fields
        del out_cell_dict, hru_type_in_dict, cell_xy_dict


        ## Flow Accumulation
        logging.info('\nCalculating initial flow accumulation')
        flow_acc_full_obj = FlowAccumulation(flow_dir_obj)
        logging.info('  Only keeping flow_acc >= {0}'.format(flow_acc_threshold))
        flow_acc_full_obj = Con(
            flow_acc_full_obj >= flow_acc_threshold, flow_acc_full_obj)
        flow_acc_full_obj.save(flow_acc_full_path)
        ##logging.info('  Only keeping active cells for subset')
        ##flow_acc_sub_obj = Con(
        ##    ((hru_type_in_obj == 1) | (hru_type_in_obj == 2)), flow_acc_full_obj)
        ##flow_acc_sub_obj.save(flow_acc_sub_path)

        ## Flow accumulation and stream link with lakes
        logging.info('\nCalculating flow accumulation & stream link (w/ lakes)')
        flow_acc_obj = Con(
            ((hru_type_in_obj == 1) | (hru_type_in_obj == 2)), flow_acc_full_obj)
        ##flow_acc_obj.save(flow_acc_sub_path)
        stream_link_obj = StreamLink(flow_acc_obj, flow_dir_obj)
        stream_link_obj.save(stream_link_a_path)
        del flow_acc_obj, stream_link_obj
        
        ## Flow accumulation and stream link without lakes
        logging.info('Calculating flow accumulation & stream link (w/o lakes)')
        flow_acc_obj = Con((hru_type_in_obj == 1), flow_acc_full_obj)
        ##flow_acc_obj.save(flow_acc_sub_path)
        stream_link_obj = StreamLink(flow_acc_obj, flow_dir_obj)
        stream_link_obj.save(stream_link_b_path)
        del flow_acc_obj, stream_link_obj

        ## Initial Stream Link
        ##logging.info('\nCalculating initial stream link')
        ##stream_link_obj = StreamLink(flow_acc_obj, flow_dir_obj)
        ##stream_link_obj.save(stream_link_path)
        ## Calculate stream link with and without lakes
        ## Initial Stream Order (w/ lakes)
        logging.info('Calculating stream order (w/ lakes)')
        logging.debug(
            '  Using SHREVE ordering so after 1st order are removed, '+
            '2nd order will only be dangles')
        stream_order_obj = StreamOrder(
            stream_link_a_path, flow_dir_obj, 'SHREVE')
        stream_order_obj.save(stream_order_path)
        ## Stream Length (cell count w/o lakes)
        logging.info('Calculating stream length (cell count w/o lakes)')
        stream_length_obj = Lookup(stream_link_b_path, 'Count')
        stream_length_obj.save(stream_length_path)
        ## Filter 1st order segments
        logging.info(
            ('\nFilter all 1st order streams with length < {0}'+
	     '\nKeep all higher order streams').format(
                 flow_length_threshold))
        ## Stream length is nodata for lakes, so put lakes back in
        ## This removes short 1st order streams off of lakes
        flow_mask_obj = (
            (hru_type_in_obj == 2) | (stream_order_obj >= 2) |
            ((stream_order_obj == 1) &
             (stream_length_obj >= flow_length_threshold)))
        flow_mask_obj.save(flow_mask_path)
        flow_acc_sub_obj = Con(flow_mask_obj, flow_acc_full_obj)
        flow_acc_sub_obj.save(flow_acc_sub_path)
        del flow_mask_obj, stream_order_obj, stream_length_obj

        ## Final Stream Link
        logging.info('\nCalculating final stream link')
        stream_link_obj = StreamLink(flow_acc_sub_obj, flow_dir_obj)
        ## Get count of streams for automatically setting lake_seg_offset
        if not lake_seg_offset:
            lake_seg_count = int(
                arcpy.GetCount_management(stream_link_obj).getOutput(0))
            n = 10 ** math.floor(math.log10(lake_seg_count))
            lake_seg_offset = int(math.ceil((lake_seg_count + 1) / n)) * int(n)
            logging.info(
                ('  lake_segment_offset was not set in the input file\n'+
                 '  Using automatic lake segment offset: {0}').format(
                     lake_seg_offset))
        elif set_lake_flag:
            logging.info(
                ('  Using manual lake segment offset: {0}').format(lake_seg_offset))
        ## Include lake cells into "stream_link" before calculating watersheds
        ## Watershed function doesn't work for negative values
        ## Convert lakes to large positive numbers for Watershed
        ## ISEG needs to be negative values though
        if set_lake_flag:
            logging.info(
                ('  Including lakes as {0} + {1}\n'+
                 '  This will allow for a watershed/subbasin for the lakes\n'+
                 '  {2} will be save as negative of {0} though').format(
                     hru.lake_id_field, lake_seg_offset, hru.iseg_field))
            stream_link_obj = Con(
                (hru_type_in_obj == 2),
                (lake_id_obj + lake_seg_offset), stream_link_obj)
        stream_link_obj.save(stream_link_path)
 
        ## Watersheds
        logging.info('Calculating watersheds')
        watersheds_obj = Watershed(flow_dir_obj, stream_link_obj)
        watersheds_obj.save(watersheds_path)
        del stream_link_obj, watersheds_obj
        
        ## Subbasins
        logging.info('Calculating subbasins')
        subbasin_obj = Watershed(
            flow_dir_obj, subbasin_points_path, subbasin_zone_field)
        subbasin_obj.save(subbasin_path)
        del subbasin_obj

        ## Basins
        logging.info('Calculating basins')
        basin_obj = Basin(flow_dir_obj)
        basin_obj.save(basin_path)
        del basin_obj


        ## Reset HRU_TYPE to 0 for inactive water cells
        logging.info('Clearing inactive water HRU_TYPE')
        hru_type_obj = Con((hru_type_in_obj == 4), 0, hru_type_in_obj)
        hru_type_obj.save(hru_type_path)
        del hru_type_obj
        del hru_type_in_obj

        ## Clear subbasin value if HRU_TYPE_IN is 0
        logging.info('Clearing subbasin ID for inactive cells')
        subbasin_obj = SetNull(
            Raster(hru_type_path), Raster(subbasin_path), "VALUE=0")
        subbasin_obj.save(subbasin_path)
        del subbasin_obj

        #### DEADBEEF
        #### Instead of clearing HRU_TYPE if no subbasin, try to keep all edge
        ####   cells in the model by giving them a common subbasin
        #### Clear HRU_TYPE for cells with no subbasin
        ##logging.info('Adjusting HRU_TYPE for cells not in subbasins')
        ##hru_type_obj = Con(
        ##    (((hru_type_in_obj == 1) | (hru_type_in_obj == 2)) &
        ##     (IsNull(subbasin_obj) == 1)), 
        ##    0, hru_type_in_obj)


        ## Stream polylines
        logging.info('Calculating stream polylines')
        ## ArcGIS fails for raster_to_x conversions on a network path
        ## You have to go through an in_memory file first
        streams_temp = os.path.join('in_memory', 'streams')
        StreamToFeature(
            stream_link_path, flow_dir_obj, streams_temp, 'NO_SIMPLIFY')
        arcpy.CopyFeatures_management(streams_temp, streams_path)
        arcpy.Delete_management(streams_temp)
        del streams_temp
        
		
        ## Write values to hru_polygon
        logging.info('\nExtracting stream parameters')
        vt_list = [
            [watersheds_path, hru.irunbound_field],
            [stream_link_path, hru.iseg_field],
            ##[flow_dir_path, hru.flow_dir_field], 
            [subbasin_path, hru.subbasin_field],
            [hru_type_path, hru.type_field]]
        mem_point_path = os.path.join('in_memory', 'hru_point')
        arcpy.CopyFeatures_management(hru.point_path, mem_point_path)
        ExtractMultiValuesToPoints(mem_point_path, vt_list, 'NONE')
        del vt_list

        ## Read values from points
        logging.info('  Reading cell values')
        data_dict = defaultdict(dict)
        fields = [
            hru.irunbound_field, hru.iseg_field, 
            hru.subbasin_field, hru.type_field, hru.fid_field]
        ##fields = [
        ##    hru.irunbound_field, hru.iseg_field, hru.flow_dir_field,
        ##    hru.subbasin_field, hru.type_field, hru.fid_field]
        with arcpy.da.SearchCursor(mem_point_path, fields) as s_cursor:
            for row in s_cursor:
                for i, field in enumerate(fields[:-1]):
                    ## Set nodata or inactive cells to 0
                    if row[i] is None or (int(row[-2]) == 0):
                        data_dict[int(row[-1])][field] = 0
                    else:
                        data_dict[int(row[-1])][field] = int(row[i])
                del row
        del fields
        
        ## ISEG for lake cells must be -1 * LAKE_ID, not LAKE_ID + OFFSET
        for k in data_dict.keys():
            irunbound = data_dict[k][hru.irunbound_field]
            iseg = data_dict[k][hru.iseg_field]
            if irunbound > lake_seg_offset:
                data_dict[k][hru.irunbound_field] = lake_seg_offset - irunbound
            if iseg > lake_seg_offset:
                data_dict[k][hru.iseg_field] = lake_seg_offset - iseg
                
        ##data_dict = dict([(k,v) for k,v in data_dict.items()])
        ## Write values to polygon
        logging.info('  Writing values to polygons')
        fields = [
            hru.irunbound_field, hru.iseg_field, 
            hru.subbasin_field, hru.type_field, hru.fid_field]
        ##fields = [
        ##    hru.irunbound_field, hru.iseg_field, hru.flow_dir_field, 
        ##    hru.subbasin_field, hru.type_field, hru.fid_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                row_dict = data_dict.get(int(row[-1]), None)
                for i, field in enumerate(fields[:-1]):
                    if row_dict:
                        row[i] = row_dict[field]
                    else:
                        row[i] = 0   
                u_cursor.updateRow(row)
                del row_dict, row
        del fields


        ## Write sink values to hru_polygon
        logging.info('\nExtracting sink values')
        mem_point_path = os.path.join('in_memory', 'hru_point')
        arcpy.CopyFeatures_management(hru.point_path, mem_point_path)
        vt_list = [
            [dem_sink8_path, hru.dem_sink8_field],
            [dem_sink4_path, hru.dem_sink4_field]]
        ExtractMultiValuesToPoints(mem_point_path, vt_list, 'NONE')
        ## Read wink values from points
        logging.info('  Reading sink values')
        data_dict = defaultdict(dict)
        fields = [hru.dem_sink8_field, hru.dem_sink4_field, hru.fid_field]
        with arcpy.da.SearchCursor(mem_point_path, fields) as s_cursor:
            for row in s_cursor:
                for i, field in enumerate(fields[:-1]):
                    ## Set nodata or inactive cells to 0
                    if row[i] is None:
                        data_dict[int(row[-1])][field] = 0
                    else:
                        data_dict[int(row[-1])][field] = float(row[i])
                del row
        ## Write sink values to polygon
        logging.info('  Writing sink values to polygons')
        fields = [hru.dem_sink8_field, hru.dem_sink4_field, hru.fid_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as u_cursor:
            for row in u_cursor:
                row_dict = data_dict.get(int(row[-1]), None)
                for i, field in enumerate(fields[:-1]):
                    if row_dict:
                        row[i] = row_dict[field]
                    else:
                        row[i] = 0   
                u_cursor.updateRow(row)
                del row_dict, row


        ## Cleanup
        arcpy.Delete_management(mem_point_path)
        del mem_point_path, vt_list, data_dict, fields


        ## Re-Calculate HRU_ELEV
        ##logging.info('Calculating HRU_ELEV from DEM_ADJ')
        ##logging.info('  Converting from meters to feet')
        ##arcpy.CalculateField_management(
        ##    hru.polygon_path, hru_elev_field,
        ##    ## Convert meters to feet
        ##    '!{0}! * 3.28084'.format(dem_adj_field), 'PYTHON')

        
        ## Cleanup
        del dem_fill_obj
        if set_lake_flag:
            del lake_id_obj
        del flow_dir_obj
        del flow_acc_full_obj
        del flow_acc_sub_obj

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

    finally:
        ##pass
        try: arcpy.CheckInExtension('Spatial')
        except: pass
        ##arcpy.ResetEnvironments()
        ##try: logging.info('\nTime elapsed = {0} seconds\n'.format(clock()-start))
        ##except: pass

################################################################################

#### Flood fill algorithm
def flood_fill(test_array, four_way_flag=True, edge_flt=None):
    input_array = np.copy(test_array)
    input_rows, input_cols = input_array.shape
    h_max = np.max(input_array * 2.0)
    logging.debug("  Hmax: %s" % (h_max/2.0))
    
    ## Since ArcGIS doesn't ship with SciPy (only numpy), don't use ndimage module
    if four_way_flag:
        el = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]).astype(np.bool)
    else:
        el = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool)
    ##if four_way_flag:
    ##    el = ndimage.generate_binary_structure(2,1).astype(np.int)
    ##else:
    ##    el = ndimage.generate_binary_structure(2,2).astype(np.int)
    
    ## Build data/inside/edge masks
    data_mask = ~np.isnan(input_array)
    
    ## Since ArcGIS doesn't ship with SciPy (only numpy), don't use ndimage module
    inside_mask = np_binary_erosion(data_mask, structure=el)
    ##inside_mask = ndimage.binary_erosion(data_mask, structure=el)

    edge_mask = (data_mask & ~inside_mask)
    ## Initialize output array as max value test_array except edges
    output_array = np.copy(input_array)
    output_array[inside_mask] = h_max
    ## Set edge pixels less than edge_flt to edge_flt
    if edge_flt:
        output_array[edge_mask & (output_array<=edge_flt)]=edge_flt

    ## Build priority queue and place edge pixels into queue
    put = heapq.heappush
    get = heapq.heappop
    fill_heap = [
        (output_array[t_row,t_col], int(t_row), int(t_col), 1)
        for t_row, t_col in np.transpose(np.where(edge_mask))]
    heapq.heapify(fill_heap)
    ##logging.info("    Queue Size: %s" % len(fill_heap))

    ## Cleanup
    del data_mask, edge_mask, el
    ##logging.info("    Prep Time: %s" % (clock()-start_total))

    while True:
        try:
            h_crt, t_row, t_col, edge_flag = get(fill_heap)
        except IndexError:
            break
        for n_row, n_col in [
            ((t_row-1),t_col), ((t_row+1),t_col),
            (t_row,(t_col-1)), (t_row,(t_col+1))]:
            ## Skip cell if outside array edges
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
    
    ## Pad output array (binary_erosion) with extra cells around the edge
    ## so that structuring element will fit without wrapping.
    ## A 3x3 structure, will need 1 additional cell around the edge
    ## A 5x5 structure, will need 2 additional cells around the edge
    output_shape = tuple(
        ss + dd - 1 for ss,dd in zip(input_array.shape, structure.shape))
    input_pad_array = np.zeros(output_shape).astype(np.bool)
    input_pad_array[1:rows+1,1:cols+1] = input_array
    binary_erosion = np.zeros(output_shape).astype(np.bool)

    ## Cast structure element to boolean
    struc_mask = structure.astype(np.bool)

    ## Iterate over each cell
    for row in xrange(rows):
        for col in xrange(cols):
            ## The value of the output pixel is the minimum value of all the
            ##   pixels in the input pixel's neighborhood.
            binary_erosion[row+1,col+1] = np.min(
                input_pad_array[row:row+3, col:col+3][struc_mask])
    return binary_erosion[1:rows+1,1:cols+1]

def raster_to_array(input_item, mask_extent=None):
    ## input_item can be a raster_obj or raster_path
    try: input_obj = Raster(input_item)
    except TypeError: input_obj = input_item
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
    ## Integer type raster can't have NaN values, will only set floats to NaN
    if (output_array.dtype == np.float32 or
        output_array.dtype == np.float64):        
        output_array[output_array == input_nodata] = np.NaN
        output_nodata = np.NaN
    else:
        output_nodata = int(input_nodata)
    return output_array, output_nodata

def array_to_raster(input_array, output_path, pnt, cs, mask_array=None):
    output_array = np.copy(input_array)
    ## Float arrays have to have nodata set to some value (-9999)
    if (output_array.dtype == np.float32 or
        output_array.dtype == np.float64):
        output_nodata = -9999
        output_array[np.isnan(output_array)] = output_nodata
    ## Boolean arrays need to be converted unsigned ints
    elif output_array.dtype == np.bool:
        output_array = output_array.astype(np.uint8)
        output_nodata = 255
    ## If a mask array is give, assume all 0 values are nodata
    if np.any(mask_array): output_array[mask_array == 0] = output_nodata
    output_obj = arcpy.NumPyArrayToRaster(
        output_array, pnt, cs, cs, output_nodata)
    output_obj.save(output_path)
    del output_obj
    arcpy.DefineProjection_management(
        output_path, env.outputCoordinateSystem)
    arcpy.CalculateStatistics_management(output_path)

################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='DEM To Streams',
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

    ## Create Basic Logger
    logging.basicConfig(level=args.loglevel, format='%(message)s')

    #### Get GSFLOW config file
    ##ini_re = re.compile('\w*.ini$', re.I)
    ##try: 
    ##    ini_path = sys.argv[1]
    ##except IndexError:
    ##    ini_path = get_ini_file(workspace, ini_re, 'gsflow_dem_parameters')
    ##del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', os.getcwd()))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    ## Convert input file to an absolute path
    if os.path.isfile(os.path.abspath(args.ini)):
        args.ini = os.path.abspath(args.ini)

    ## Calculate GSFLOW Flow Parameters
    gsflow_flow_parameters(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
