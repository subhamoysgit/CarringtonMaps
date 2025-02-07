import matplotlib.pyplot as plt
import numpy as np
from sunpy.map import Map
import astropy.units as u
from astropy.io import fits
import scipy.ndimage
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import wcs_to_celestial_frame


def hmi_hg_proj(f,resolution,angle):
    """Heliographic projection

    Args:
        f (str): filename
        resolution (float): deg/pixel
        angle (float): heliocentric angle

    Returns:
        HMImap: sunpy.map
        hg (np.ndarray): heliographic map
        hg_strk (np.ndarray): heliographic streak map
    """    
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
    wcs1 = wcs_to_celestial_frame(HMImap.wcs)
    hg = np.ones((180*res,180*res))*(-1000)
    hg_strk = hg.copy()
    n_lat = np.linspace(-90,89,180*res)
    n_lon = np.linspace(-90,89,180*res)
    lon,lat = np.meshgrid(n_lon,n_lat)
    sc = SkyCoord(lon*u.deg, lat*u.deg, frame='heliographic_stonyhurst',obstime = HMImap.date,observer=wcs1.observer)
    pix_x = (HMImap.meta['crpix1'] + sc.transform_to("helioprojective").Tx/(HMImap.meta['cdelt1']*u.arcsec)).astype('int')
    pix_y = (HMImap.meta['crpix2'] + sc.transform_to("helioprojective").Ty/(HMImap.meta['cdelt2']*u.arcsec)).astype('int')
    for i in range(180*res):
        for j in range(180*res):
            hg[i,j] = conv[pix_y[i,j],pix_x[i,j]]
            phi = (np.pi/180.)*(j/res)
            hg_strk[i, j] = (hg[i, j]>-1)*np.cos(phi - np.pi/2)**4
    return HMImap,hg, hg_strk


def hmi_car_proj(f,resolution):
    """Carrington projection

    Args:
        f (str): filename
        resolution (float): deg/pixel

    Returns:
        carr (np.ndarray): Carrington map
        carr_strk (np.ndarray): Carrington streak map
    """    
    HMImap, hg, hg_strk = hmi_hg_proj(f,resolution,80)
    res = int(1/resolution)
    crlon_obs = HMImap.meta['crln_obs']
    npix = int((crlon_obs-180)/(15/90))
    carr = -1000*np.ones((180*res,360*res))
    carr_strk = carr.copy()
    carr[:,90*res:270*res] = hg
    carr_strk[:,90*res:270*res] = hg_strk
    carr = scipy.ndimage.shift(carr, np.array([0,npix]),mode='wrap')
    carr_strk = scipy.ndimage.shift(carr_strk, np.array([0, npix]),mode='wrap')
    return carr, carr_strk