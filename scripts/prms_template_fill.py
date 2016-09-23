#--------------------------------
# Name:         prms_template_fill.py
# Purpose:      Fill PRMS Parameter File Template
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2016-02-26
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
# import csv
import datetime as dt
# import itertools
import logging
import operator
import os
# import re
import sys

import arcpy
# from arcpy import env

from support_functions import *


def prms_template_fill(config_path, overwrite_flag=False, debug_flag=False):
    """Fill PRMS Parameter Template File

    Args:
        config_file (str): Project config file path
        ovewrite_flag (bool): if True, overwrite existing files
        debug_flag (bool): if True, enable debug level logging

    Returns:
        None
    """

    # Initialize hru_parameters class
    hru = HRUParameters(config_path)

    # Open input parameter config file
    config = ConfigParser.ConfigParser()
    try:
        config.readfp(open(config_path))
    except:
        logging.error(('\nERROR: Config file could not be read, ' +
                       'is not an input file, or does not exist\n' +
                       'ERROR: config_file = {0}\n').format(config_path))
        sys.exit()

    # Log DEBUG to file
    log_file_name = 'prms_template_log.txt'
    log_console = logging.FileHandler(
        filename=os.path.join(hru.log_ws, log_file_name), mode='w')
    log_console.setLevel(logging.DEBUG)
    log_console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(log_console)
    logging.info('\nFilling PRMS Parameter File Template')

    # Read parameters from config file
    hru.fid_field = config.get('INPUTS', 'orig_fid_field')
    prms_parameter_path = config.get('INPUTS', 'prms_parameter_path')
    prms_dimen_csv_path = config.get('INPUTS', 'prms_dimen_csv_path')
    prms_param_csv_path = config.get('INPUTS', 'prms_param_csv_path')
    parameter_ws = config.get('INPUTS', 'parameter_folder')

    # Scratch workspace
    try:
        scratch_name = config.get('INPUTS', 'scratch_name')
    except:
        scratch_name = 'in_memory'

    # Strings to search PRMS parameter file for
    file_header_str = 'Default file generated by model\nVersion: 1.7'
    dimen_header_str = '** Dimensions **'
    param_header_str = '** Parameters **'
    break_str = '####'

    # Check input paths
    if not arcpy.Exists(hru.polygon_path):
        logging.error(
            '\nERROR: HRU Shapefile ({0}) does not exist'.format(
                hru.polygon_path))
        sys.exit()
         
    if not os.path.isfile(prms_dimen_csv_path):
        logging.error('\nERROR: The dimensions CSV file does not exist\n')
        sys.exit()
    
    if not os.path.isfile(prms_param_csv_path):
        logging.error('\nERROR: The parameters CSV file does not exist\n')
        sys.exit()
    
    if os.path.isfile(prms_parameter_path):
        os.remove(prms_parameter_path)

    #Get the total number of HRUs from shapefile
    hru_count = int(arcpy.GetCount_management(hru.polygon_path).getOutput(0))

    # Read in Parameter Dimensions CSV file
    logging.info('\nReading PRMS dimensions CSV file')
    dimen_size_dict = dict()
    with open(prms_dimen_csv_path, 'r') as input_f:
        dimen_lines = input_f.readlines()
    input_f.close()
    
    # Dimensions can be set to a value, a field, or not set
    dimen_lines = [l.strip().split(',') for l in dimen_lines]
    header = dimen_lines[0]
    for line in dimen_lines[1:]:
        dimen_size = line[header.index('SIZE')]
        if dimen_size =='CALCULATED':
            pass
        elif not dimen_size:
            dimen_size_dict[line[header.index('NAME')]] = ''
        else:
            dimen_size_dict[line[header.index('NAME')]] = int(dimen_size)
        del dimen_size
        
    # These parameters equal the total number of HRUs  in the HRU shapefile
    for dimen_name in ['ngw','nhru', 'nhrucell', 'nssr']:
        dimen_size_dict[dimen_name] = hru_count
        logging.info('  {0} = {1}'.format(
            dimen_name, dimen_size_dict[dimen_name]))

    # Getting number of lakes
#     logging.info('\nCalculating number of lake cells')
#     logging.info('  Lake cells are {0} >= 0'.format(
#         hru.lake_id_field))
#     value_fields = (hru.id_field, hru.lake_id_field)
#     with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
#         dimen_size_dict['nlake'] = len(list(
#             [int(row[1]) for row in s_cursor if int(row[1]) > 0]))
#     logging.info('  nlakes = {0}'.format(dimen_size_dict['nlake']))
#
#     # Getting number of stream cells
#     logging.info('Calculating number of stream cells')
#     logging.info('  Stream cells are {0} >= 0'.format(
#         hru.krch_field))
#     value_fields = (hru.id_field, hru.krch_field)
#     with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
#         dimen_size_dict['nreach'] = len(list(
#             [int(row[1]) for row in s_cursor if int(row[1]) > 0]))
#     logging.info('  nreach = {0}'.format(dimen_size_dict['nreach']))
#    
    # Getting number of stream segments
    logging.info('Calculating number of unique stream segments')
    stream_segments = arcpy.da.UpdateCursor(hru.stream_path, ["OBJECTID"])
    nsegment = 0
    for segment in stream_segments:
        nsegment +=1
    dimen_size_dict['nsegment']=nsegment
    logging.info('  nsegment = {0}'.format(dimen_size_dict['nsegment']))

    # Getting number of subbasins
    logging.info('Calculating number of unique subbasins')
    logging.info('  Subbasins are {0} >= 0'.format(
        hru.subbasin_field))
    value_fields = (hru.id_field, hru.subbasin_field)
    with arcpy.da.SearchCursor(hru.polygon_path, value_fields) as s_cursor:
        dimen_size_dict['nsub'] = len(list(set(
            [int(row[1]) for row in s_cursor if int(row[1]) > 0])))
    logging.info('  nsub = {0}'.format(dimen_size_dict['nsub']))

    # Link HRU field names to parameter names in '.param'
    param_name_dict = dict()
    param_width_dict = dict()
    param_dimen_count_dict = dict()
    param_dimen_names_dict = dict()
    param_values_count_dict = dict()
    param_type_dict = dict()
    param_default_dict = dict()
    param_values_dict = defaultdict(dict)

    # Read in parameters from CSV
    logging.info('\nReading PRMS parameters CSV file')
    with open(prms_param_csv_path, 'r') as input_f:
        param_lines = input_f.readlines()
    input_f.close()
    param_lines = [l.strip().split(',') for l in param_lines]
    header = param_lines[0]
    for line in param_lines[1:]:
        # Get parameters from CSV line
        param_name = line[header.index('NAME')]
        print param_name
        param_width = line[header.index('WIDTH')]
        # This assumes multiple dimensions are separated by semicolon
        dimen_names = line[header.index('DIMENSION_NAMES')].split(';')

        # Check that parameter type is 1, 2, 3, or 4
        param_type = int(line[header.index('TYPE')])
        if param_type not in [1, 2, 3, 4]:
            logging.error(
                ('\nERROR: Parameter type {0} is invalid' +
                 '\nERROR: {1}').format(param_type, line))
            sys.exit()
        # This will initially read defaults in as a list
        param_default = line[header.index('DEFAULT_VALUE'):]
        # Removing empty strings avoids checking ints/floats
        param_default = [val for val in param_default if val]
        # For empty lists, set to none
        if not param_default:
            param_default = None
        # For single value lists, get first value
        # Check that param_default is a number or field name
        elif len(param_default) == 1:
            param_default = param_default[0]
            #Check for current type =float when expecting integer
            if isfloat(param_default) and param_type == 1:
                param_default = int(param_default)
            #Check for current type =float when expecting float or double
            elif isfloat(param_default) and param_type in [2, 3]:
                param_default = float(param_default)
            elif param_default == 'CALCULATED':
                pass

            elif arcpy.ListFields(hru.polygon_path, param_default):
                pass
            else:
                logging.error(
                    ('\nERROR: Default value {0} was not parsed' +
                     '\nERROR: {1}').format(param_default, line))
                sys.exit()

        # For multi-value lists, convert values to int/float
        elif len(param_default) >= 2:
            #Expecting integers
            if param_type == 1:
                param_default = map(int, param_default)
            #Expecting float/double
            elif param_type in [2, 3]:
                param_default = map(float, param_default)
            else:
                logging.error(
                    ('\nERROR: Default value {0} was not parsed' +
                     '\nERROR: {1}').format(param_default, line))
                sys.exit()

        # Check that dimension names are valid
        for dimension in dimen_names:
            if dimension not in dimen_size_dict.keys():
                logging.error(
                    ('\nERROR: The dimension {0} is not set in the ' +
                     'dimension CSV file').format(dimension))
                sys.exit()

        # Calculate number of dimensions
        dimen_count = str(len(dimen_names))

        # Calculate number of values
        values_count = prod( [int(dimen_size_dict[dn]) for dn in dimen_names if dimen_size_dict[dn]])

        # Write parameter to dictionaries
        param_name_dict[param_name] = param_name  #What is the point of dict where keys = values? Seems unecessary..
        param_width_dict[param_name] = param_width
        param_dimen_count_dict[param_name] = dimen_count
        param_dimen_names_dict[param_name] = dimen_names
        param_values_count_dict[param_name] = values_count
        param_type_dict[param_name] = param_type
        param_default_dict[param_name] = param_default
    #END of PRMS Parameter and dimensions CSV read in

    # Apply default values to full dimension of default parameters
    logging.info('\nSetting static parameters from defaults')
    for param_name, param_default in param_default_dict.items():
        param_values_count = param_values_count_dict[param_name]
        # Skip if not set
        if param_default is None:
            continue
        # Skip if still a string (field names)
        elif type(param_default) is str:
            continue
        # For float/int, apply default across dimension size
        elif type(param_default) is float or type(param_default) is int:
            for i in range(param_values_count):
                param_values_dict[param_name][i] = param_default
        # For lists of floats, match up one-to-one for now
        elif len(param_default) == param_values_count:
            for i in xrange(param_values_count):
                param_values_dict[param_name][i] = param_default[i]
        else:
            logging.error('\nERROR: The default value(s) ({0}) could not be ' +
                 'broadcast to the dimension length ({1})').format(
                     param_default, param_values_count)
            sys.exit()
    
    # Begin to read in HRU parameter data from Hru shapefile
    logging.info('\nReading in variable parameters from HRU shapefile')
   
    #Take all the parameters that are defined using string value to another value and create a dict
    param_field_dict = {}
    hru_param_field_dict = {}
    strm_param_field_dict = {}

    for key,value in param_default_dict.items():
        
        #Stream parameters have to come from stream shapefile not hru, so collect the related params here and index later.
        if type(value) is str and value in ["TOSEGMENT","X_COEF","K_COEF","OBSIN_SEGMENT"]:
            strm_param_field_dict[key]=value
        # Add all string valued params except "calculated" 
        elif type(value) is str and value != "CALCULATED":
            param_field_dict[key] = value
    
     
    arc_value_fields = param_field_dict.values()
    strm_arc_value_fields = strm_param_field_dict.values()

    # Use an ID to uniquely identify each cell as to place its value correctly under the key
    identifier = hru.id_field
    if  identifier not in arc_value_fields:
        arc_value_fields.append(identifier)
    
    # Read in each cell parameter value
    s_cursor = arcpy.da.SearchCursor(hru.polygon_path, arc_value_fields)
    for row in s_cursor:
        #Use the identifier to uniquely assign each value in cursor
        param_row_id = row[arc_value_fields.index(identifier)]

        #Iterate through and add all values in the row to our dict
        for param_name,arc_param_name in param_field_dict.items():
            param_values_dict[param_name][param_row_id] = row[arc_value_fields.index(arc_param_name)]
            #logging.debug('{0}, {1}, {2}'.format(row[arc_value_fields.index(identifier)],param_name,row[arc_value_fields.index(arc_param_name)]))
    del s_cursor
    

    identifier = "OBJECTID"
    if  identifier not in arc_value_fields:
        strm_arc_value_fields.append(identifier)
    #Now add the stream parameters from the stream shape file.
    stream_cursor = arcpy.da.SearchCursor(hru.stream_path, strm_arc_value_fields)
    for row in stream_cursor:
        #Use the identifier to uniquely assign each value in cursor
        param_row_id = row[strm_arc_value_fields.index(identifier)]

        #Iterate through and add all values in the row to our dict
        for param_name,arc_param_name in strm_param_field_dict.items():
            param_values_dict[param_name][param_row_id] = row[strm_arc_value_fields.index(arc_param_name)]
    del stream_cursor
    
    # Calculate mean monthly maximum temperature for all active cells
    logging.info('\nCalculating tmax_index')
    logging.info('  Converting Celsius to Farenheit')
    param_name_dict['tmax_index'] = 'tmax_index'
    param_width_dict['tmax_index'] = 15
    param_dimen_count_dict['tmax_index'] = 1
    param_dimen_names_dict['tmax_index'] = ['nmonths']
    param_values_count_dict['tmax_index'] = dimen_size_dict['nmonths']
    param_type_dict['tmax_index'] = 2
    tmax_field_list = ['TMAX_{0:02d}'.format(m) for m in range(1, 13)]
    for i, tmax_field in enumerate(tmax_field_list):
        tmax_values = [row[1] for row in arcpy.da.SearchCursor(
            hru.polygon_path, (hru.type_field, tmax_field),
            where_clause='"{0}" >= 1'.format(hru.type_field))]
        tmax_c = sum(tmax_values) / len(tmax_values)
        tmax_f = 1.8 * tmax_c + 32
        param_values_dict['tmax_index'][i] = tmax_f
        logging.info('  {0} = {1}'.format(
            tmax_field, param_values_dict['tmax_index'][i]))
        del tmax_values

    cell_dict = dict()
    fields = [
        hru.type_field, hru.krch_field, hru.lake_id_field,
        hru.subbasin_field, hru.flow_dir_field]
#         ,hru.col_field, hru.row_field, hru.id_field]
    for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
        # Skip inactive cells
        if int(row[0]) == 0:
            continue
        # Skip non-lake and non-stream cells
        elif (int(row[1]) == 0 and int(row[2]) == 0):
            continue
        # Read in parameters
        cell = (int(row[5]), int(row[6]))
        # next_row_col(FLOW_DIR, CELL)
        # HRU_ID, SUBBASIN, NEXT_CELL
        cell_dict[cell] = [
            int(row[7]), int(row[3]), next_row_col(int(row[4]), cell)]
        del cell

    # # DEADBEEF - lake_hru is not used in PRMS 3.0.X or gsflow
    # #   It is used in PRMS 4.0 though
    # # lake_hru parameter
    # logging.info('\nCalculating LAKE_HRU from HRU_ID for all lake HRU\'s')
    # param_name_dict['lake_hru'] = 'lake_hru'
    # param_width_dict['lake_hru'] = 0
    # param_dimen_count_dict['lake_hru'] = 1
    # param_dimen_names_dict['lake_hru'] = ['nlake']
    # param_values_count_dict['lake_hru'] = dimen_size_dict['nlake']
    # param_type_dict['lake_hru'] = 1
    # lake_hru_id_list = [
    #    row[1] for row in arcpy.da.SearchCursor(
    #        hru.polygon_path, (hru.type_field, hru.id_field))
    #    if int(row[0]) == 2]
    # for i,lake_hru_id in enumerate(sorted(lake_hru_id_list)):
    #    # logging.debug('  {0} {1}'.format(i, lake_hru_id))
    #    param_values_dict['lake_hru'][i] = lake_hru_id


    # # Add lake HRU's to groundwater cascades
    # logging.info('Modifying CRT groundwater parameters for all lake HRU\'s')
    # logging.info('  gw_up_id = HRU_ID (lake)')
    # logging.info('  gw_down_id = 0')
    # # logging.info('  gw_strmseg_down_id = OUTSEG')
    # logging.info('  gw_strmseg_down_id = 2')
    # logging.info('  gw_pct_up = 1')
    # field_list = [hru.type_field, hru.id_field, hru.outseg_field,
    #              hru.outflow_field]
    # lake_hru_id_dict = dict([
    #    (row[1], row[2])
    #    for row in arcpy.da.SearchCursor(hru.polygon_path, field_list)
    #    if int(row[0]) == 2 and int(row[3]) == 0])
    # for lake_hru_id, outseg in sorted(lake_hru_id_dict.items()):
    #    if lake_hru_id == 9128:
    #        print lake_hru_id, outseg
    #    # raw_input('ENTER')
    #    i = dimen_size_dict['ncascdgw']
    #    dimen_size_dict['ncascdgw'] += 1
    #    param_values_dict['gw_up_id'][i] = lake_hru_id
    #    param_values_dict['gw_down_id'][i] = 0
    #    # DEADBEEF - PRMS didn't like when set to OUTSEG, but 2 worked?
    #    param_values_dict['gw_strmseg_down_id'][i] = outseg
    #    # param_values_dict['gw_strmseg_down_id'][i] = 2
    #    # DEADBEEF - Trying 0
    #    param_values_dict['gw_strmseg_down_id'][i] = 0
    #    param_values_dict['gw_pct_up'][i] = 1.00
    #    # print param_values_dict['gw_up_id'][i]
    #    # print param_values_dict['gw_down_id'][i]
    #    # print param_values_dict['gw_strmseg_down_id'][i]
    #    # print param_values_dict['gw_pct_up'][i]
    # param_values_count_dict['gw_up_id'] = dimen_size_dict['ncascdgw']
    # param_values_count_dict['gw_down_id'] = dimen_size_dict['ncascdgw']
    # param_values_count_dict['gw_strmseg_down_id'] = dimen_size_dict['ncascdgw']
    # param_values_count_dict['gw_pct_up'] = dimen_size_dict['ncascdgw']
    # logging.info('  ncascade = {0}'.format(dimen_size_dict['ncascade']))
    # logging.info('  ncascdgw = {0}'.format(dimen_size_dict['ncascdgw']))
    # raw_input('ENTER')
 
    # Write dimensions/parameters to PRMS param file
    logging.info('\nWriting parameter file')
    with open(prms_parameter_path, 'w') as output_f:
        output_f.write(file_header_str + '\n')
        # Dimensions
        output_f.write(dimen_header_str + '\n')
        # Write dimensions that are known first
        # remove_list = []
        logging.info('  Set dimensions')
        for dimen_name, dimen_size in sorted(dimen_size_dict.items()):
            if not dimen_size:
                continue
            logging.debug('    {0}'.format(dimen_name))
            output_f.write(break_str+'\n')
            output_f.write(dimen_name+'\n')
            output_f.write(str(dimen_size)+'\n')
            # DEADBEEF - It seems bad to remove items during iteration
            del dimen_size_dict[dimen_name]
            # remove_list.append(dimen_name)
        # for dimen_name in remove_list: del dimen_size_dict[dimen_name]
        # Then write unset dimensions
        logging.info('  Unset dimensions')
        for dimen_name in sorted(dimen_size_dict.keys()):
            logging.debug('  {0}'.format(dimen_name))
            output_f.write(break_str+'\n')
            output_f.write(dimen_name+'\n')
            output_f.write(str(dimen_size_dict[dimen_name])+'\n')

        # Parameters
        output_f.write(param_header_str + '\n')
        # Write unset parameters first
        logging.info('  Unset parameters')
        for param_name in sorted(param_name_dict.keys()):
            if param_name in param_values_dict.keys():
                continue
            logging.debug('    {0}'.format(param_name))
            output_f.write(break_str+'\n')
            output_f.write('{0} {1}\n'.format(
                param_name, param_width_dict[param_name]))
            output_f.write('{0}\n'.format(param_dimen_count_dict[param_name]))
            for dimen_name in param_dimen_names_dict[param_name]:
                output_f.write(dimen_name + '\n')
            output_f.write(str(param_values_count_dict[param_name]) + '\n')
            param_type = param_type_dict[param_name]
            output_f.write(str(param_type) + '\n')
            output_f.write('' + '\n')
            # DEADBEEF - It seems bad to remove items during iteration
            del param_name_dict[param_name]
            # del param_width_dict[param_name]
            # del param_dimen_count_dict[param_name]
            # del param_dimen_names_dict[param_name]
            # del param_values_count_dict[param_name]
            # del param_type_dict[param_name]

        # Then write set parameters
        logging.info('  Set parameters')
        for param_name in sorted(param_name_dict.keys()):
            logging.debug('  {0}'.format(param_name))
            output_f.write(break_str+'\n')
            output_f.write('{0} {1}\n'.format(
                param_name, param_width_dict[param_name]))
            output_f.write('{0}\n'.format(param_dimen_count_dict[param_name]))
            for dimen_name in param_dimen_names_dict[param_name]:
                output_f.write(dimen_name + '\n')
            output_f.write(str(param_values_count_dict[param_name]) + '\n')
            param_type = param_type_dict[param_name]
            output_f.write(str(param_type) + '\n')
            for i, param_value in param_values_dict[param_name].items():
                if param_type == 1:
                    output_f.write('{0:d}'.format(int(param_value)) + '\n')
                elif param_type == 2:
                    output_f.write('{0:f}'.format(param_value) + '\n')
                elif param_type == 3:
                    output_f.write('{0:f}'.format(param_value) + '\n')
                elif param_type == 4:
                    output_f.write('{0}'.format(param_value) + '\n')

    # Close file
    output_f.close()
    logging.info('\nDone!')

def prod(iterable):
    return reduce(operator.mul, iterable, 1)


def isfloat(s):
    """"""
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False

# class dimension():
#    def __init__(self, i, data_lines):
#        self.i = i
#        self.NAME = data_lines[0]
#        self.SIZE = int(data_lines[1])
#        print self.NAME, self.SIZE

# class parameter():
#    # type_dict = dict()
#    # type_dict[1] = 'INTEGER'
#    # type_dict[2] = 'FLOAT'
#    # type_dict[3] = 'DOUBLE'
#    # type_dict[4] = 'STRING'
#    def __init__(self, i, data_lines):
#        self.i = i
#        # # Not all names have a width (hvr_hru_pct, gvr_hru_id, gvr_cell_id)
#        try: self.NAME, self.WIDTH = data_lines[0].split()
#        except ValueError: self.NAME, self.WIDTH = data_lines[0], 0
#        # # There can be multiple dimensions
#        self.NO_DIMENSIONS = int(data_lines[1])
#        self.DIMENSION_NAMES = []
#        for i in range(self.NO_DIMENSIONS):
#            self.DIMENSION_NAMES.append(data_lines[2+i])
#        self.N_VALUES = int(data_lines[3+i])
#        self.TYPE = data_lines[4+i]
#        self.VALUE = data_lines[5+i:]
#        print self.NAME, self.WIDTH, self.NO_DIMENSIONS,
#        print self.DIMENSION_NAMES, self.N_VALUES, self.TYPE


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

    # Fill PRMS Parameter Template File
    prms_template_fill(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)
