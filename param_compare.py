#--------------------------------
# Name:         param_compare.py
# Purpose:      compare the outputted para file with another one from a different group
# Author:       Micah Johnson
# Created       2016-10-20
# Python:       2.7
#--------------------------------


class ParamFile(object):
    def __init__(self,filename,**kwargs):
        if type(filename)==str:
            self.filename = filename
        else:
            ValueError("Expected filename type str not {0} ".format(type(filename)))
            sys.exit()
        self.param_dict = {}
        params = {}
        
        with open(self.filename,'r') as f:
            param_lines = f.readlines()
        
        #Go through and find every parameter name
        while i < len(param_lines):
            line = param_lines[i]
            
            if "####" in line:
                print param_lines[i+1]
            

if __name__ == '__main__':
    our_param = ParamFile('C:\Users\Public\Documents\Micahs_RCEW\prms_input\25k_rcew.param')