# -*- coding: utf-8 -*-
"""
Created on Thu Jul 09 11:16:42 2015

@author: smudd
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os

# This function takes a string and looks for the seperators
# it then reformats to the correct operating system
def ReformatSeperators(path):
    
    path = RemoveEscapeCharacters(path)
    #print "Path is: " + path    
    
    # loook for the various seperators in the string
    if "\\" in path:
        #print "I found \\"
        path = RemoveEscapeCharacters(path)
        splitpath = path.split("\\")
    elif "/" in path:
        #print "I found ''/''"
        if "//" in path:
            #print "I found ''//''"
            splitpath = path.split("//")
            
        else:
            splitpath = path.split("/")
    else:
        print("I did not find a valid seperator")
        splitpath = "NULL"

    # Now reconstitute the path using the operating system seperator
    # get the size of the split path
    n_spaces = len(splitpath)   
    newpath = ""   
    for s in range (0,n_spaces-1):
        newpath = newpath+splitpath[s]+os.sep 
            
    # now the last element
    newpath = newpath+splitpath[n_spaces-1]
    
    newpath2 = RemoveEscapeCharacters(newpath)
    
    return newpath2    

# This function takes a path with any seperators and converts the seperators
# to the current operating system and adds a seperator at the end of the
# path                    
def AppendSepToDirectoryPath(path):
    
    #first check to see if there are escape characters
    RemoveEscapeCharacters(path)    
    
    # now reformat the seperators
    newpath = ReformatSeperators(path)
    
    # now check to see if the last character is a sep
    if newpath[-1] != os.sep:
        newpath = newpath+os.sep
    return newpath

# This function takes a filename (that could include a full directory path)
# and lops off the path as well as the extension
def GetFilePrefix(filename):
    newfilename = ReformatSeperators(filename)
    
    # split the file
    splitname = newfilename.split(os.sep)
    fname = splitname[-1]
    splitfname = fname.split(".")
    
    fileprefix = splitfname[0]
    return fileprefix

# This gets the filename without the path
def GetFileNameNoPath(filename):
    newfilename = ReformatSeperators(filename)
    
    # split the file
    splitname = newfilename.split(os.sep)
    fname = splitname[-1]
    return fname


# This gets the last directory level
# say if the path is home\yo\ma\
# then this function returns ma
def GetLastDirectoryLevel(path):
    newpathname = ReformatSeperators(path)
    pathname = AppendSepToDirectoryPath(newpathname)
    
    # now split the path
    splitpath = pathname.split(os.sep)
    
    # now remove the final path level
    n_spaces = len(splitpath)   
    newpath = splitpath[n_spaces-2] 

    return newpath
        
        
# This gets the path of a filename (with a full directory path)
def GetPath(filename):
    newfilename = ReformatSeperators(filename)
    
    splitname = newfilename.split(os.sep)
    n_levels = len(splitname)   
    newpath = ""   
    for s in range (0,n_levels-1):
        newpath = newpath+splitname[s]+os.sep
        
    return newpath

# This function finds the level of the path
def GetPathLevel(path):
    newpathname = ReformatSeperators(path)
    pathname = AppendSepToDirectoryPath(newpathname)
    splitname = pathname.split(os.sep)
    n_levels = len(splitname)  
    return n_levels
    
    
# This removes one level from the directory    
def RemoveDirectoryLevel(path):
    # Make sure the formatting is appropriate and there is a seperator at the end
    newpathname = ReformatSeperators(path)
    pathname = AppendSepToDirectoryPath(newpathname)
    
    # now split the path
    splitpath = pathname.split(os.sep)
    
    # now remove the final path level
    n_spaces = len(splitpath)   
    newpath = "" 
    for s in range (0,n_spaces-2):
        newpath = newpath+splitpath[s]+os.sep
        
    return newpath
    


# This is necessary because of stupid windows seperators    
def RemoveEscapeCharacters(line):
    line = line.rstrip().replace('\n', '\\n')
    line = line.rstrip().replace('\b', '\\b') 
    line = line.rstrip().replace('\r', '\\r')  
    line = line.rstrip().replace('\t', '\\t')  
    line = line.rstrip().replace('\n', '\\n')  
    line = line.rstrip().replace('\f', '\\f')  
    line = line.rstrip().replace('\v', '\\v') 
        
    # this last one deals with the infuriating special case of \b
    line = line.rstrip().replace('\x08', '\\b')     
    
    # get rid of leading and trailing spaces
    line = line.strip()     
    
    return line
    
def RemoveWhitespace(line):
    line = line.replace(" ", "")
    return line

# This function takes a string. If it is an integer, it returns an integer, 
# if it is floating point, it returns that, otherwise it returns the string.     
def ParseStringToType(A_string):
    try:
        return int(A_string)
    except ValueError:
        try: 
            return float(A_string)
        except ValueError:
            return A_string    

# This takes a list and returns either ints, floats or strings
# If any element in the list is a string, all the elements are converted to string
# If there are no strings, it is converted to ins unless a single element is a
# float, in which case all elements are floats        
def ParseListToType(A_list):   
    have_string = False
    have_float = False
    new_list = []    
    
    for element in A_list:
        try:
            int(element)
        
        except ValueError:
            
            try: 
                float(element)
                have_float = True
            
            except ValueError:
                have_string = True
                
    if have_string:
        #print "I found a string"
        for element in A_list:
            new_list.append(element)
    else:
        if have_float:
            #print "These are floats"
            for element in A_list:
                new_list.append(float(element))
        else:
            #print "These are ints"
            for element in A_list:
                new_list.append(int(element)) 
                
    return new_list
                
    
        
        
