"""
This is a tester function for some of the mapping scripts

Author: SMM

Date 05/01/2018
"""





import sys
import os



def set_path():
    """
    This ensures that the path is included to the LSDMappingTools packages.
    It avoids a complicated pip installation but you need the LSDTT environment working. 
    
    Author: SMM
    Date 08/01/2018
    
    """
    # We need to import the path since we don't install mapping tools directly
    this_dir = os.getcwd()
    print("This dir is: "+this_dir)

    LSDMT_dir = this_dir.split(os.sep)

    LSDMT_dir = LSDMT_dir[:-1]
    LSDMT_path = (os.sep).join(LSDMT_dir)

    print("The path is: "+LSDMT_path)
    sys.path.append(LSDMT_path)
    

def print_welcome():
    """
    This displays a welcome screen.
    
    Author: SMM and FJC
    Date: 08/01/2018
    """

    print("\n\n=======================================================================")
    print("Hello! I'm going to test if your LSDTopoTools python enviroment is working.")
    print("You need to tell me what analysis to do:")
    print("1 = Basic hillshade map")
    print("2 = Map with some basins")
    print("3 = Map wit a channel network")
    print("Run with: python Mapping_test.py 1 ")
    print("The outputs can be found as files in this folder.")
    print("=======================================================================\n\n ")

    

def main(argv):
    """
    The main driving function.
    
    Author: SMM
    Date: 08/01/2018
    """
    
    if len(argv) == 0:
        print_welcome()
        quit()
    else:
        print("Let me load the LSDMappingTools functions for you.")
        set_path()
    
 
    print("Let me check the fonts. They are:")
    import matplotlib.font_manager as fm
    flist = fm.get_fontconfig_fonts()
    names = [fm.FontProperties(fname=fname).get_family() for fname in flist]   
    print(names)
    
    fm.findfont('Liberation Sans',rebuild_if_missing=True)

    print("The arguments are: ")
    print(argv)
    
    # import Maping tools
    import LSDMapWrappers as LSDMW
    DataDir = os.getcwd()+os.sep
    DataFname = "WA"
    
    argv[0] = int(argv[0])

    if argv[0] == 1:
        print("Getting basic hillshade")
        LSDMW.SimpleHillshade(DataDir,DataFname)
    elif argv[0] == 2:
        print("Plotting some basins")
        LSDMW.PrintBasins(DataDir,DataFname)
    elif argv[0] == 3:
        print("Plotting the channels")
        LSDMW.PrintChannels(DataDir,DataFname) 
    else:
        print("I didn't understand what you wanted.")
        print("Your choice was:" + str(argv[0]))
    



#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])