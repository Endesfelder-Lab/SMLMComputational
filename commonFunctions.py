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
    
def create_trajectories(tracked_data,frame1,frame2,maxDistance,tracksCounter):
  #Make sub-matrix of this frame and the next frame
  framematrix = tracked_data[tracked_data[:,0]==frame1,:];
  nextframematrix = tracked_data[tracked_data[:,0]==frame2,:];

  #Now we check that there are localizations in both of these frames - if one
  #of the frames does not have localizations, we cannot try to track the
  #localizations. Note that we already know that in this case, frame 0 and
  #frame 1 have localizations, but this will not always be the case.
  if (len(framematrix) > 0 and len(nextframematrix) > 0):
    #Now we can find all the nearest neighbours between all localizations
    # on this frame and on the next frame
    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(nextframematrix[:,1:3])
    NearestNeighbors(n_neighbors=1)
    foundnn = neigh.kneighbors(framematrix[:,1:3])
    foundnn = np.asarray(foundnn)

    # Now we have to select only those entries that are < maxDistance. We want
    # to extract the IDs of the localizations in frame n+1 that belong to them.
    NeighbourIDs = foundnn[:,foundnn[1] < maxDistance][1].astype(int)
    # We also want to extract the ID of the original localizations by
    # finding the row index of the nearest neighbours.
    OriginIDs = np.where(foundnn[1] < maxDistance)[0].astype(int)

    # For every found neighbour, we make both neighbours the same track-id.
    # First, we check that the neighbour on frame n has or doesn't have an
    # track-id yet, then we set the neighbour on frame n+1
    # to the same value, or to a new track-id if none is assigned yet
    # We loop over all found IDs
    for i in range(0,len(NeighbourIDs)):
      #We get the localization-ID of the neighbour in frame n+1
      neighbourID = NeighbourIDs[i]
      #We also get the localization-ID of the original localization in frame n
      originID = OriginIDs[i]
      #Prevent linkage if the neighbour in frame n+1 is already linked -
      #this will not happen in this example, but it might happen when
      #skipping frames in later modules.
      if nextframematrix[neighbourID,4] == 0:
        # We check that the localization that will be included in a trajectory
        # is not yet part of an existing trajectory - if it is part of
        # an existing track, the index in column 4 will be higher than 0
        if framematrix[originID,4] == 0:
          #If it's not linked yet, set it and the neighbour to a new track-id value
          framematrix[originID,4] = tracksCounter
          nextframematrix[neighbourID,4] = tracksCounter
          tracksCounter += 1
        else:
          #If it is linked, set it to the track-id value of the neighbour in frame n
          nextframematrix[neighbourID,4] = framematrix[originID,4]

        #Finally, we  provide the distance to the next emitter in the
        #track on the next frame
        framematrix[originID][5] = foundnn[0][originID][0]

  #Now we have fully filled frame matrix and nextframematrix variables, but
  #we need to fill these back in into the original tracked_data matrix. We do
  #this by looking up the values via the idCounter value
  if len(framematrix) > 0:
      tracked_data[tracked_data[:,0]==frame1] = framematrix
  if len(nextframematrix) > 0:
      tracked_data[tracked_data[:,0]==frame2] = nextframematrix
      
  #Return the required parameters
  return tracked_data, tracksCounter

