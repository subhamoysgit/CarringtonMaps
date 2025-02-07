#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 17:27:27 2021

@author: schatterjee
"""

import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import cv2
import sunpy.map
from sunpy.map import Map
import astropy.units as u
import numpy as np
import os
from astropy.io import fits
import os
import drms

import requests
from bs4 import BeautifulSoup
import julian
from datetime import datetime
import numpy as np
import scipy.ndimage

import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.coordinates

from astropy.wcs.utils import wcs_to_celestial_frame, custom_wcs_to_frame_mappings
from sunpy.coordinates import Helioprojective

def hmi_hg_proj(f,resolution,angle):
  res = int(1/resolution)
  HMI_fits = fits.open(f)
  HMI_fits.verify('fix')
  HMImap = Map(HMI_fits[1].data, HMI_fits[1].header)
  HMImap = HMImap.rotate(recenter=True)
  x, y = np.meshgrid(*[np.arange(v.value) for v in HMImap.dimensions]) * u.pixel
  hpc_coords = HMImap.pixel_to_world(x, y)
  rSun = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / HMImap.rsun_obs
  conv = HMImap.data/1.3
  conv[rSun>np.sin(np.pi*angle/180)] = -1000
  conv =  np.clip(conv,-1000,1000)
  HMImap = Map(conv,HMImap.meta)
  wcs = HMImap.wcs
  wcs1 = wcs_to_celestial_frame(HMImap.wcs)
  hg = np.ones((180*res,180*res))*(-1000)
  n_lat = np.linspace(-90,89,180*res)
  n_lon = np.linspace(-90,89,180*res)
  lon,lat = np.meshgrid(n_lon,n_lat)
  sc = SkyCoord(lon*u.deg, lat*u.deg, frame='heliographic_stonyhurst',obstime = HMImap.date,observer=wcs1.observer)
  pix_x = (HMImap.meta['crpix1'] + sc.transform_to("helioprojective").Tx/(HMImap.meta['cdelt1']*u.arcsec)).astype('int')
  pix_y = (HMImap.meta['crpix2'] + sc.transform_to("helioprojective").Ty/(HMImap.meta['cdelt2']*u.arcsec)).astype('int')
  for i in range(180*res):
    for j in range(180*res):
      hg[i,j] = conv[pix_y[i,j],pix_x[i,j]]
  return HMImap,hg


def hmi_car_proj(f,resolution):
  HMImap,hg = hmi_hg_proj(f,resolution,80)
  res = int(1/resolution)
  crlon_obs = HMImap.meta['crln_obs']
  npix = int((crlon_obs-180)/(15/90))
  carr = -1000*np.ones((180*res,360*res))
  carr[:,90*res:270*res] = hg
  carr = scipy.ndimage.shift(carr, np.array([0,npix]),mode='wrap')
  return carr

import pickle
import os
from datetime import datetime
import julian
year = input('enter year : ')
year = str(year)
months = {'2010':'12','2011':'11','2012':'10','2013':'09','2014':'08','2015':'07','2016':'06','2017':'05','2018':'04','2019':'03','2020':'02','2021':'01'}
str_0 = year + months[year] + '01_0000'
#print(str_0)
dt_0 = datetime(int(str_0[:4]),int(str_0[4:6]),int(str_0[6:8]),int(str_0[9:11]),int(str_0[11:13]),0,0)
jd_0 = julian.to_jd(dt_0,fmt='jd')
time = []
st = []
carr_map = np.zeros((1080,2160,500))
k = 0
for root,dirs,files in os.walk('/home/schatterjee/Documents/'+year):
  for name in files:
    str_n = name[11:24] 
    st.append(str_n)
    dt = datetime(int(str_n[:4]),int(str_n[4:6]),int(str_n[6:8]),int(str_n[9:11]),int(str_n[11:13]),0,0)
    jd = julian.to_jd(dt, fmt='jd')
    ex = int(np.round(24*3600*(jd-jd_0)))
    time.append(ex)
    carr_map[:,:,k] = hmi_car_proj('/home/schatterjee/Documents/'+year+'/'+name,(15/90))
    k = k+1
    print(str_n)
    print(k)

p = [carr_map[:,:,:k],st,ex]
pickle.dump(p,open('/home/schatterjee/Documents/'+year+'/carr_map_video.p','wb'),protocol=4)
