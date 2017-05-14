"""This class contains calculation for sorting data in order to make the plotting scripts lighter
@author Boris Gailleton"""
import numpy as np

def source_num_of_longest_river(arr):
    """
    This function will return the source number of the longuest channel by counting the number of occurence of each sources number.
    It may be an issue if your channel have similar size, I will create a fancier function that calculate the actual lenght of each channel if required.
    """
    counount = 0
    max_count = -999
    max_count_index = -999

    for i in range(1,arr.shape[0]):
        if(arr[i] == arr[i-1]):
            counount += 1
        else:
            if (counount > max_count):
                max_count_index = arr[i-1]
                max_count = counount
            counount = 0

    return max_count_index
