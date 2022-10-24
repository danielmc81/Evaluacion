#!/usr/bin/env python
# coding: utf-8

import cartopy.io.img_tiles as cimgt
from wrf import to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,cartopy_ylim, latlon_coords, ALL_TIMES
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import matplotlib
from netCDF4 import Dataset
import cartopy.crs as ccrs
import numpy as np
from scipy.io import netcdf_file
from datetime import datetime, timedelta
import pandas as pd
import sys
import multiprocessing as mp
import glob
import os

def llminmax(var):
    lats, lons = latlon_coords(var)

    lat_min = to_np(lats.min()) + 0.0
    lat_max = to_np(lats.max()) - 0.0
    lon_min = to_np(lons.min()) + 0.0
    lon_max = to_np(lons.max()) - 0.0
    
    extent = [lon_max, lon_min, lat_max, lat_min]
    
    return lats, lons, extent

def mapa(lats, lons, extent, var, request, cmap, levels, title, ftitle, currtime):
    sx, sy = 14, 14  
    dpi = 72
    alpha = 0.3
    
    plt.figure(figsize=(sx, sy))
    ax = plt.axes(projection=request.crs)

    t_hmo = pd.to_datetime(currtime) - timedelta(hours=7)
    t_utc = pd.to_datetime(currtime)

    ax.set_extent(extent, crs=ccrs.PlateCarree())

    norm = matplotlib.colors.BoundaryNorm(levels, cmap.N)

    cf = plt.contourf(to_np(lons), to_np(lats), to_np(var), 
                 transform=ccrs.PlateCarree(), norm=norm,
                 cmap=cmap,levels=levels, alpha=alpha)
    
    if dominio == "d01":
        shrink=0.92
        zoom=6      
    if dominio == "d02":
        shrink= 0.99
        zoom=8
    if dominio == "d03":
        shrink= 0.99
        zoom=8        
    if dominio == "d04":
        shrink= 0.80
        zoom=8              
    if dominio == "d05":
        shrink=1
        zoom=7
    if dominio == "d06":
        shrink=1
        zoom=9
    if dominio == "d07":
        shrink=0.86
        zoom=8
    if dominio == "d08":
        shrink=0.91
        zoom=8
    
    plt.colorbar(ax=ax, shrink=shrink, orientation='horizontal', pad=0.025)
    
    ax.set_title(title, fontsize=18)
    ax.set_title("Tiempo de Hermosillo\n" + str(t_hmo), loc='left', fontsize=14)
    ax.set_title("UTC\n" + str(t_utc), loc='right', fontsize=14)
    
    filename = ftitle + "-" + str(currtime) + "_" + dominio + ".png"
    ax.add_image(request, zoom)
    plt.show()
    plt.savefig(filename)
    plt.close()

def humedad_relativa(lats, lons, extent, rh, request, time):
    cmap = get_cmap("jet")
    levels = np.arange(0, 101, 1)
    title = "Humedad Relativa (%)"
    ftitle = "humedad_relativa"

    mapa(lats, lons, extent, rh, request, cmap, levels, title, ftitle, time)

def temperatura_2m(lats, lons, extent, t2, request, time):
    cmap = get_cmap("hsv").reversed()
    levels = np.arange(0, 50, 1)
    title = r"Temperatura 2m ($\degree$C)"
    ftitle = "temp2m"
    
    mapa(lats, lons, extent, t2, request, cmap, levels, title, ftitle, time)

def agua_precipitable(lats, lons, extent, pw, request, time):
    cmap = get_cmap("terrain").reversed()
    levels = np.arange(0, 80, 1)
    title = "Vapor de Agua Precipitable\n (mm)"
    ftitle = "agua_precip"
    
    mapa(lats, lons, extent, pw, request, cmap, levels, title, ftitle, time)

def precipitacion(lats, lons, extent, p, title, time):
    levels = [0.1, 0.5, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
    precip_colors = [
        "#04e9e7", 
        "#019ff4",  
        "#0300f4",  
        "#02fd02",  
        "#01c501",  
        "#008e00",  
        "#fdf802",  
        "#e5bc00",  
        "#fd9500",  
        "#fd0000",  
        "#d40000",  
        "#bc0000",  
        "#f800fd",  
        "#9854c6",  
        "#fdfdfd"]  

    precip_colormap = matplotlib.colors.ListedColormap(precip_colors)
    cmap = precip_colormap

    title = "Precipitación 1hr (mm)"
    ftitle = "precip"
    
    mapa(lats, lons, extent, p, request, cmap, levels, title, ftitle, time)   
    
def precipitacion_acumulada(lats, lons, extent, pa, request, time):
    levels = [0.1, 0.5, 1, 2, 5, 10, 15, 20, 30, 40, 50, 80, 100, 200, 300]
    precip_colors = [
        "#04e9e7", 
        "#019ff4",
        "#0300f4",
        "#02fd02",
        "#01c501",
        "#008e00", 
        "#fdf802",  
        "#e5bc00", 
        "#fd9500",
        "#fd0000", 
        "#d40000", 
        "#bc0000",
        "#f800fd",
        "#9854c6",
        "#fdfdfd"]  

    precip_colormap = matplotlib.colors.ListedColormap(precip_colors)
    cmap = precip_colormap
    
    title = "Precipitación Acumulada (mm)"
    ftitle = "precip_acum"
    
    mapa(lats, lons, extent, pa, request, cmap, levels, title, ftitle, time)

def radiacion_solar(lats, lons, extent, rs, request, time):
    cmap = get_cmap("hot").reversed()
    levels = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]    

    title = "Radiación Solar (w/m2)"
    ftitle = "rad_solar"
    
    mapa(lats, lons, extent, rs, request, cmap, levels, title, ftitle, time)


def plot():
    print("Plotting Temperatura")
    for time in range(len(times)):
        p = mp.Process(target=temperatura_2m, args=(lats, lons, extent, t2[time][0], request, times[time].values))        
        p.start()     
    p.join() 

    print("Plotting Precipitacion Acumulada")
    for time in range(len(times)):
        p = mp.Process(target=precipitacion_acumulada, args=(lats, lons, extent, (rainc[time]+rainnc[time]), request, times[time].values))        
        p.start()        
    p.join()

    print("Plotting Precipitacion")
    for time in range(len(times)):  
        if time > 0:      
            p = mp.Process(target=precipitacion, args=(lats, lons, extent, ((rainc[time]+rainnc[time])-(rainc[time-1]+rainnc[time-1])), request, times[time].values))        
        else:
            p = mp.Process(target=precipitacion, args=(lats, lons, extent, (rainc[time]+rainnc[time]), request, times[time].values))        
        p.start()       
    p.join()

    print("Plotting Vapor de Agua Precipitable")
    for time in range(len(times)):
        p = mp.Process(target=agua_precipitable, args=(lats, lons, extent, pw[time], request, times[time].values))        
        p.start()        
    p.join()

    print("Radiacion Solar")
    for time in range(len(times)):
        p = mp.Process(target=radiacion_solar, args=(lats, lons, extent, rs[time], request, times[time].values))        
        p.start()        
    p.join()

    print("Humedad Relativa")
    for time in range(len(times)):
        p = mp.Process(target=humedad_relativa, args=(lats, lons, extent, rh[time], request, times[time].values))        
        p.start() 
    p.join()

if __name__ == '__main__':
    ruta =  sys.argv[1] + "/wrfout*"

    request = cimgt.GoogleTiles()

    for wrfout in glob.glob(ruta):
        url = wrfout
        fname = os.path.basename(wrfout)
        dominio = fname[7:10]
                
        ds = Dataset(url)
        
        print("Inicializando ", dominio)
        times = getvar(ds, "times", timeidx=ALL_TIMES)        
        rh = getvar(ds, "rh2", timeidx=ALL_TIMES)
        t2 = getvar(ds, "tc", timeidx=ALL_TIMES)
        pw = getvar(ds, "pw", timeidx=ALL_TIMES)
        rainc = ds["RAINC"]
        rainnc = ds["RAINNC"]
        rs = ds["SWDOWN"]

        lats, lons, extent = llminmax(t2)

        plot()  