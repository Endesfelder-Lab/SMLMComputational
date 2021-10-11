# Main import module
import numpy as np # Used in modules 1, 2, 3, 4a, 4b, 5, 6, 7, 8, 9a, 9b
import pandas as pd # Used in modules 4a, 4b, 5, 6, 7, 8, 9a, 9b
import matplotlib.pyplot as plt # Used in modules 1, 2, 4a, 4b, 5, 6, 7, 8, 9a, 9b  
from sklearn.neighbors import NearestNeighbors # Used in modules 3, 4a, 7, 8, 9b
import sys # Used in module 4a, 8
from IPython.display import clear_output # Used in modules 4a, 4b, 6, 9b
from scipy.interpolate import interp1d # Used in module 4b
import cv2 # Used in modules 2, 5
from scipy.io import loadmat # Used in modules 5
from skimage import io # Used in modules 1, 2
from scipy.signal import convolve2d # Used in module 6
from scipy.spatial import distance # Used in module 3, 7
from scipy.optimize import curve_fit # Used in modules 7, 9b
from sklearn.cluster import DBSCAN # Used in module 8
from scipy.ndimage import gaussian_filter # Used in module 2
from skimage.feature.peak import peak_local_max # Used in module 2

# Functions to read datafiles:
# Output will always be an array with columns frame - x - y - intensity
# Input can be ThunderSTORM file (readCSV) or rapidSTORM file (readTXT)
def readCSV(filename):
  import pandas
  #Determine which fields we want to read
  fields = ['x [nm]', 'y [nm]','frame', 'intensity [photon]']

  #Read the csv with only specific headers
  data = pandas.read_csv(filename, usecols=['frame','x [nm]','y [nm]','intensity [photon]'])
  #Convert the data to a regular array
  data = data.to_numpy()
  return data

def readTXT(filename):
  import pandas
  #Read the TXT as a CSV file, with tab-separated data, using only columns 1-4
  data = pandas.read_csv(filename, header=None, skiprows=1, sep=" ", usecols=[0,1,2,3])
  #Convert the data to a regular array
  data = data.to_numpy()
  #Reorder to frame-x-y-int for consistency
  data = data[:,[2,0,1,3]]
  return data

def progress_bar(iterated_object, progressbar_len=50):
    progress = (i+1)/len(iterated_object)
    block = int(progressbar_len*progress)
    clear_output(wait=True)
    text = "Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (progressbar_len - block), progress * 100)
    print(text)
