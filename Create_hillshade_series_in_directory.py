"""
This function reads hillshades and elevation files in a directory


INSTRUCTIONS FOR ffmpeg
This will produce a series of images that you can convert to a movie using ffmpeg

Here is a typical command line
ffmpeg -framerate 5 -pattern_type glob -i '*.png' -vcodec libx264 -s 1230x566 -pix_fmt yuv420p movie7.mp4

Some gotchas:
The -framerate flag must come after ffmpeg (I think)
You can be more specific about the pattern but this seems to work
For the -s flag you must have pixel sizes that are divisible by 2.

Written by Simon Mudd
June 2017
git
GPL3
"""


import LSDChiMappingExamples as CME
import LSDMapWrappers as MW
import os
from sys import platform


def get_filenames(root):
    # Create and empty list for the filenames

    print("The root is: "+root)

    these_filenames = []
    n = 0
    while 1:
        filename = "%s%d.bil" % (root, n+1)
        print("The filename is: "+filename)
        if os.path.isfile(filename):
             n += 1
             these_filenames.append(filename)
        else:
            return these_filenames

def run_plots(DataDirectory,Base_file):

    root = DataDirectory+Base_file
    filenames = get_filenames(root)
    counter = 0
    for fname in filenames:

        counter = counter+1
        # remove the extension form the files
        this_base_file = fname[:-4]
        if "win" in platform:
            split_file = this_base_file.split("\\")
        elif "linux" in platform:
            split_file = this_base_file.split("/")
        print("Split file is: ")
        print split_file
        this_base_file = split_file[-1]
        print("This base file is: "+ this_base_file)


        MW.SimpleHillshadeForAnimation(DataDirectory,this_base_file,cmap = "terrain", dpi = 250, imgnumber=counter, full_basefile = root)



if __name__ == "__main__":

    """
    This is just a few lines for keeping track of how long the program is taking.
    You can ignore it.
    """
    import time
    tic = time.clock()

    """
    These lines tell the example functions where your files are. If you are using
    the recommended file setup you won't need to modify these lines.

    If you have set up your own directory structure you will need to modify
    the directory names.
    """
    #DataDirectory = "S:\\movern_analysis\hokkaido\\"
    #Base_file = "hokkaido_points"
    #DataDirectory = "T:\\analysis_for_papers\\Xian\\"
    #DataDirectory = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\LSDTT_chi_examples\\'
    #Base_file = 'Xian2'
    #DataDirectory = 'S:\\movern_analysis\\model_runs\\muddpile\\big_fluvial_test\\'
    DataDirectory = "/home/s0923330/muddpile_test/test1/"
    Base_file = 'LSDRM'

    """
    These lines are used to run the examples. Each line calls a function in the directory
    ../LSDChiMappingExamples/

    The individual examples are different python scripts which you can inspect at your
    lesiure. By playing with these scripts you can learn how to use our plotting tools.

    To run the examples simply comment or uncomment these lines by adding or
    removing the comment symbol, #, below.
    """

    run_plots(DataDirectory,Base_file)

    toc = time.clock()
    print("This took: "+str(toc - tic)+" units of time")
