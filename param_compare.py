
from fileinput import filename
import __main__

def get_local_params(param_lines,f_len):
    """
    Search a PRMS params file that was generated locally.
     
    *Param Indicator = ####
     
    *File is a single column
     
    """
    i=0
    names = []
    while i < f_len:
        line = param_lines[i-1]
        #search for all the delimiters.
        if "####" in line:
            name = param_lines[i].split(" ")[0]
            names.append(name.split('\n')[0])
        i+=1
    return names

def get_interface_params(param_lines,f_len):
    """
    Search a PRMS params file that was generated locally.
     
    *Param Indicator = @
     
     
    *data is multidimensional
    """
    i=1
    names = [] 
    while i < f_len:
        line = param_lines[i]
        #search for all the delimiters.
        if "@" in line:
            lst_name = param_lines[i].split(",")[1:-1]
            for param in lst_name:
                #Have we selected this param yet?
                if param not in names:
                    try:
                        #found param values not strings of names
                        float(param)
                    except:
                        #New param name, grab it.
                        names.append(param)
            
        i+=1
    return names
    
class ParamFile(object):
    def __init__(self,filename,**kwargs):
        if type(filename) == str:
            self.filename = filename
        else:
            ValueError("\nExpected filename to be type str not {0}".format(type(filename)))
        
        
        #Try opening the file.
        with open(self.filename,'r') as f:
            param_lines = f.readlines()
            f_len = len(param_lines)
        i=1
        #Is the param file from the online interface or generated locally?

        #First line @s,parameters
        if param_lines[0][0]=="@":
            print "\nWeb interface generated file detected."
            self.names  = get_interface_params(param_lines,f_len)

        #First line text from user.
        else:
            print "\nLocally generated file detected."
            self.names  = get_local_params(param_lines,f_len)

        
#             if i>0:
#                 if param_lines[i-1][0:3]=="####":
#                     print param_lines[i]
#  
def print_missing_names(their_lst, our_lst):
    for their_name in their_lst:
        for our_name in our_lst:
            if their_name==our_name:
                break
        else:
            print their_name
if __name__=='__main__':
    my_file = "C:/Users/Public/Documents/Micahs_RCEW/prms_input/25k_rcew.params"
    their_file = "C:/Users/Public/Documents/George_online_files/params.csv"
    #their_file = "C:/Users/Public/Documents/PRMS_FromBruce_scott/Tuol.scott.params"
    our_params = ParamFile(my_file)
    their_params = ParamFile(their_file)
    print "\nThese are the parameters not found in the file comparison:\n\n\t{0} \n\t\t\t\tVS \n\t{1}".format(my_file,their_file)   
    print_missing_names(their_params.names,our_params.names)        
   
