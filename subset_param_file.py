"""
Subset a single parameter file or put multiple files back together

20170407 Scott Havens 
"""


import argparse
import sys
import numpy as np
from select_param_editor import ParamEdit
import re


class Subset(ParamEdit):
    
    def __init__(self, file_name, upper=None, lower=None, selector=None, subbasin=None):
        """
        Extension of ParamEdit but can be more general
        """
        
#         if upper:
#             upper = float(upper)
#         if lower:
#             lower = float(lower)  
        if subbasin:
            subbasin = [int(s) for s in subbasin.split(',')]
            selector = 'hru_subbasin'
            parameter = 'hru_subbasin'       
                    
        self.selector = selector
        self.file = file_name
        self.parameter = parameter
        self.upper = upper
        self.lower = lower
        self.subbasin = subbasin
        self.new_value = None
        
        self.selected_col = 'all'
       
        print("\nOpening Parameter file for editing...")

        parameter_floc, ids, param_lines = self.return_edit_info()
        
        d = self.filter_by_id(param_lines, ids)
        
                        
        print("\n Done!")
        
        
    def filter_by_id(self, param_lines, hru_id):
        """
        filter the parameter file by the hru_id
        """
        
        print("{} HRU's were matched".format(len(hru_id)))
        param_lines
        new_lines = []
        
        pat = re.compile('nhru|ngw|nssr|nsub')
        it = iter(param_lines)
        dim = False
        param = False
        dims = {}
        
        
        # first get the dimensions
        for dimension_end,line in enumerate(param_lines):
            if "** Parameters **" in line:
                break
        
        i = 0
        while i < dimension_end:
            line = next(it)
            if "####" in line:
                dim_name = next(it).strip()
                dim_val = int(next(it).strip())
                dims[dim_name] = dim_val
                i += 2
            i += 1
            
        dims['nhru'] = len(hru_id)
        dims['ngw'] = len(hru_id)
        dims['nssr'] = len(hru_id)
        dims['nsub'] = len(self.subbasin)
             
        # now go look through the parameters
        it = iter(param_lines)
        while True:
#         for i,line in enumerate(param_lines):
            
            line = next(it)
            
            # check where we are in the parameter file
            if 'Dimensions' in line:
                dim = True
            if 'Parameters' in line:
                param = True
                dim = False
             
            
            # check if it's a new paramter and gather some prelim info
            if param:
                if '####' in line:
                    new_lines.append(line)
                    new_lines.append(next(it)) # parameter name
                    
                    ndims = next(it)
                    new_lines.append(ndims)
                    
                    p = 1
                    for i in range(int(ndims.strip())):
                        d = next(it)
                        new_lines.append(d)
                        
                        p *= dims[d.strip()]
                    
                    new_lines.append("{}\n".format(p))
                    
                        
                        
            
            # search for lines that have nhru, ngw, nssr, nsub
            if pat.match(line):
                
                
                new_lines.append(line) # the dimension variable
                
                if 'nsub' in line:
                    new_line = '{}\n'.format(len(self.subbasin))
                else:
                    new_line = '{}\n'.format(len(hru_id))
                    
                old_idx = int(next(it).strip())  # the old value
                new_lines.append(new_line) # the new value
                                   
                if param:

                    if 'nsub' in line:
                        idx = self.subbasin
                    else:
                        idx = hru_id
                    
                    for i in range(old_idx):
                        line = next(it)
                        if i in idx:
                            new_lines.append(line)
                    
#                     new_lines.append(line)
                    
                
            else:
                new_lines.append(line)
        
        return new_lines


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
                        help = "The parameter to be edited")
    
    parser.add_argument('--param_file','-f',
                        type=str,
                        required=True,
                        help = "The parameter file to be edited")
    

    parser.add_argument('--subbasin','-b',
                type=str,
                help = "subbasin allows a user to specify changes to specific subbasin using changes for a parameter by any other of the edits.")
    
    args = parser.parse_args()
    
  
    edit = Subset(file_name = args.param_file,
                     selector = args.by_variable,
                     upper = args.upper_thresh,
                     lower = args.lower_thresh,
                     subbasin = args.subbasin)