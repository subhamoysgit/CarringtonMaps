import pickle
import os
from datetime import datetime
import julian
import matplotlib.pyplot as plt
import numpy as np
from sunpy.map import Map
import astropy.units as u
from astropy.io import fits
import scipy.ndimage
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import wcs_to_celestial_frame
from utils import hmi_car_proj

def main():
  DIR = '/home/schatterjee/Documents/'
  year = input('enter year : ')
  year = str(year)
  months = {'2010':'12','2016':'06','2021':'01'}
  str_0 = year + months[year] + '01_0000'
  dt_0 = datetime(int(str_0[:4]),int(str_0[4:6]),
                  int(str_0[6:8]),int(str_0[9:11]),
                  int(str_0[11:13]),0,0)
  jd_0 = julian.to_jd(dt_0,fmt='jd')
  time = []
  st = []
  carr_map = np.zeros((1080,2160,500))
  carr_map_strk = carr_map.copy()
  k = 0
  resolution = (15/90)
  for _, _, files in os.walk(DIR + year):
    for name in files:
      str_n = name[11:24] 
      st.append(str_n)
      dt = datetime(int(str_n[:4]),
                    int(str_n[4:6]),
                    int(str_n[6:8]),
                    int(str_n[9:11]),
                    int(str_n[11:13]),0,0)
      jd = julian.to_jd(dt, fmt='jd')
      ex = int(np.round(24*3600*(jd-jd_0)))
      time.append(ex)
      carr, carr_strk = hmi_car_proj(DIR+year+'/'+name, resolution)
      carr_map[:,:,k] = carr
      carr_map_strk[:, :, k] = carr_strk
      k = k+1

  p = [carr_map[:,:,:k],carr_map_strk[:,:,:k],st,ex]
  pickle.dump(p,open(DIR+year+'/carr_map_video.p','wb'),protocol=4)
  # Sum carr_map*carr_strk and carr_strk over a Carrinton Rotation (CR) period, 
  # then divide the former with latter to produce the final map for a CR
  
if __name__ == '__main__':
  main()
