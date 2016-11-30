'''
2016 Micah Johnson

Allows a user to easily edit a PRMS parameter file using thresholds and specifying most values.
'''

import pandas as pd
import os, csv
import argparse
import sys


class ParamEdit(object):
    def __init__(self,file, parameter, new_value, upper=None,lower = None, selector = None ):
        
        
        if selector ==None:
            selector = parameter             
        print "\nOpening Parameter file for editing..."
        param_lines, lines_len = self.get_lines(file)
        if upper:
            upper = float(upper)
        if lower:
            lower = float(lower)    
        parameter_floc = self.find_parameter(parameter,param_lines, lines_len)
        selector_floc = self.find_parameter(selector,param_lines, lines_len)
        
        ids  = self.get_ids(selector_floc, param_lines, upper,lower)
        
        new_lines = self.set_param_values(parameter_floc, ids, new_value, param_lines)
        
        self.write_lines(file, new_lines)
        print "\n Done!"
        
    def get_lines(self,filename):
        """
        Gets all the lines from param file and returns them as a list 
        of strings and the lenght of the list.
        """
        
        #Open parameter file
        with open(filename,'r+') as f:
            lines = f.readlines()
            f.close()
        return lines, len(lines)
    
    def write_lines(self,filename,lines):
        """
        Writes all lines to file.
        """
        with open(filename,'w+') as f:
            f.writelines(lines)
            f.close()
        print "\nParameter was written to: \n{0}".format(filename)
        
            
            
    def find_parameter(self, name, lines, lines_len):
        """
        Finds a parameter and returns a given line number in the file.
        """
        i = 0
        line = lines[0]
        #Search for the parameter name.
        while name not in line and i < lines_len-1:
            i+=1
            line = lines[i]

        #EOF?
        if i ==lines_len-1:
            print "\nError: Parameter {0} was not found in parameter file.".format(name)
            self.recomend_str(name, lines)
            sys.exit()
        else:
            return i
        
    def get_ids(self, name_floc, lines, upper = None, lower = None):
        """
        Return the hru ids of hrus whose specified parameter falls between upper and lower
        
        arguments:
            name_floc       The file line location of a parameter that is being used to select hru IDS for editing
            uppper          The upper bound on the parameter being used for selecting.
            lower           The lower threshold for which the selection parameter is measured.

        returns:
            List of integers representing hru IDs.
        """
        print "\nRetrieving IDs from {0}...".format(lines[name_floc].strip())
        dimension, rows, cols, data_len, data_type, data_start, data_end = self.get_param_info(name_floc,lines)
        ids = []
        message = "\tCollecting all hru ids using parameter {0}".format(self.string_to_data(str, lines[name_floc],lines)) 
        if lower:
            message += " >= {0}".format(lower)
            if upper:
                message+= " and"
        if upper:
            message += " <= {0}".format(upper)

        print message
        
        for i in range(data_start, data_start + rows):
            value = self.string_to_data(data_type,lines[i],lines)
            #Hru id is one based.
            hru_id = i - data_start+1
            #Two cases, 1. Upper and Lower defined, 2. Upper defined lower is not.
            if upper:
                if lower:
                    if value <= upper and value >= lower:
                        ids.append(hru_id)
                else:
                    if value <= upper:
                        ids.append(hru_id)
            #Alternate case is upper is not defined and lower is.
            else:
                if lower != None:
                    if value >= lower:
                        ids.append(hru_id)
                else:
                    ids.append(hru_id)
        return ids
    
    def string_to_data(self,data_type,str_data,lines):
        """
        Convert our string data to actual data given the structure of a params file.
        """
        data = str_data.strip()
   
        if data_type == int or data_type == 1:
            try:
                return int(data)
            except:
                #Looks for a dimension instead
                return self.get_dimension(str_data,lines)
        
        elif data_type ==float or data_type in [2,3]:
            return float(data)
        
        elif data_type == str or data_type ==4:
            return data
        else:
            ValueError("Attempted to convert value {0} of type {1} to string".format(str_data,data_type))
            

    def get_dimension(self,dimension_str,lines):
        """
        Retrieves/returns dimension values who are named by a string
        """

        for line in lines:
            if "** Parameters **" in line:
                dimension_end = lines.index(line)

        print "\tSearching for dimension {0}'s value in parameter file...".format(dimension_str.strip())
        result = None
        
        for i in range(0,dimension_end+1):
            line = lines[i]
            if dimension_str in line:
                #Retrvieve the value underneath lhe line found
                result = self.string_to_data(1, lines[i+1], lines)
                break
        if result:
            return result
        else:
            print "\nERROR: Unable to find dimension {0}! Check parameter file.".format(dimension_str.strip())
            sys.exit()
                        
         
    def get_param_info(self,floc,lines):
        """
        Param file is structure such that there is the param name 
        #followed by columns, dimension of each column, total number of values, and number types
        before values actually start.
        """
        info = []
        
        #Get parameter dimensions. if the dim is more than 1 then the length of dimension data is differtent.
        dimension = self.string_to_data(1,lines[floc+1],lines)

        if dimension <2:
            info_len = 4
            for i in range(1,info_len+1):
                data = self.string_to_data(1,lines[floc+i],lines)
                info.append(data)
            cols = 1
            data_len = info[2]
            data_type = info[3]
            data_start = floc + 4
            data_end = floc+data_len
             
        else:
            info_len = 5
            for i in range(1,info_len+1):
                data = self.string_to_data(1,lines[floc+i],lines)
                info.append(data)
            cols = info[2]
            data_len = info[3]
            data_type = info[4]
            data_start = floc + 5
            data_end = floc+data_len
        
        #Rows is always the same line
        rows = info[1]
        
        if rows*cols != data_len:
            ValueError("Something went wrong when retrieving parameter information for {0}".format(lines[floc].strip()))
        else:        
            return dimension, rows, cols, data_len, data_type, data_start, data_end
    
    def recomend_str(self, search, lines):
        """
        Searchs for a 80% match of your inputted str and recommends and that match that close.
        This is meant to be used when parameter names are hard to remember.
        """
        recommendations = []
        str = self.string_to_data(4, search,lines)
        for line in lines:
            i=0
            new_str = self.string_to_data(4, line,lines)
            altered_new_str = new_str
            for s in str:
                if s in altered_new_str:
                    i+=1
                    #Removed the character if match is found
                    altered_new_str.replace(s,"",1)
            if float(i)/len(str) >0.8:
                recommendations.append(new_str)
        
        print "\nHINT:You are searching for {0}, did you mean to use these?".format(str)
        print "\n"
        for rec in recommendations:
            print "{0}".format(rec)
            
    def set_param_values(self, name_floc, hru_ids, new_value, lines):
        """
        Goes to a parameter and set the parameter values to new value according to the hru_id. returns the new list of lines edited.
        """
        print "\nSetting new values for {0}...".format(lines[name_floc].strip())         
    
        dimension, rows, cols, data_len, data_type, data_start, data_end = self.get_param_info(name_floc,lines)
        #Cycle through every columns and every line of data.
        iter =0
        prms_types = [[int,1], [float,2,3],[str,4]]
        
#         #Check user data type entry
#         for single_type in prms_types:
#             if str(type(new_value)) in single_type:
#                 break
#             
#         print "\nERROR: New value of type {0} does nto match parameter {1} which is of type {2}".format(str(type(new_value)),lines[name_floc].strip(), data_type)
#         sys.exit()
        for col in range(cols):
            #Data is stored in a single column but is specified as multiple columns, meaning the data will by cols X dimensions long.
            for i in hru_ids:
                line_id = i + data_start
                lines[line_id] = str(new_value) + '\n'
                iter +=1
            data_start+=rows

        
        print "\n\t{0} values in {1} have been changed to {2}".format(iter, lines[name_floc].strip(), new_value)       
        if iter != len(hru_ids) * cols:
            print "\n WARNING: Something appears to have gone wrong, {0} values were requested to be changed but {1} were.".format(cols*len(hru_ids), iter)
        return lines 

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= "Allows users to easily edit a parameter via the command line")
    
    parser.add_argument('--upper_thresh','-u', 
                        type = str, 
                        help = "The upper bounding value of the variable a parameter is edited by. If not set the threshold is unbounded in the upper direction")
    
    parser.add_argument('--lower_thresh','-l',
                         type = str, 
                         help = "The lower bounding value of the variable a parameter is edited by. If not set, the threshold is unbounded in the lower direction")
    
    parser.add_argument('--by_variable', '-v',
                        type = str,
                        help = "The variable we want to select parameters to edit by")

    parser.add_argument('--edit_parameter','-p',
                        type=str,
                        required=True,
                        help = "The parameter to be edited")
    
    parser.add_argument('--param_file','-f',
                        type=str,
                        required=True,
                        help = "The parameter file to be edited")
    
    parser.add_argument('--new_value','-n',
                        type=float,
                        required=True,
                        help = "The new value to be entered into all the values selected")
    
    args = parser.parse_args()
    
  
    edit = ParamEdit(file = args.param_file,
                     parameter = args.edit_parameter,
                     selector = args.by_variable,
                     upper = args.upper_thresh,
                     lower = args.lower_thresh,
                     new_value = args.new_value)
      
