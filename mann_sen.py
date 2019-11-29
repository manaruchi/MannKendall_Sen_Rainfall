# -*- coding: utf-8 -*-
"""
Application of Mann_Kendall and Sen Slope to IMD rainfall raster

Author : Manaruchi
"""
from osgeo import gdal, osr
import numpy as np
import scipy.stats as st
import sys
from tqdm import tqdm
#-------------------Input------------------------------------------------------
inpfol = r"D:\IMD\test"
outfol = r"E:\Manaruchi\python-Projects\MannKendall\out"
confidence = 0.95
#------------------------------------------------------------------------------

def mann_kendall(vals, confidence = 0.95):
    n = len(vals)
    
    box = np.ones((len(vals), len(vals)))
    box = box * 5
    sumval = 0
    for r in range(len(vals)):
        for c in range(len(vals)):
            if(r > c):
                if(vals[r] > vals[c]):
                    box[r,c] = 1
                    sumval = sumval + 1
                elif(vals[r] < vals[c]):
                    box[r,c] = -1
                    sumval = sumval - 1
                else:
                    box[r,c] = 0
    
    freq = 0
    #Lets caluclate frequency now
    tp = np.unique(vals, return_counts = True)
    for tpx in range(len(tp[0])):
        if(tp[1][tpx]>1):
            tp1 = tp[1][tpx]
            sev = tp1 * (tp1 - 1) * (2 * tp1 + 5)
            freq = freq + sev
        
    se = ((n * (n-1) * (2 * n + 5) - freq) / 18.0) ** 0.5
    
    #Lets calc the z value
    if(sumval > 0):
        z = (sumval - 1) / se
    else:
        z = (sumval + 1) / se
    
    #lets see the p value
    
    p = 2 * st.norm.cdf(-abs(z))
    

    
    #trend type
    if(p<(1-confidence) and z < 0):
        tr_type = -1
    elif(p<(1-confidence) and z > 0):
        tr_type = +1
    else:
        tr_type = 0
    
    return z, p, tr_type

def sen_slope(vals, confidence = 0.95):
    alpha = 1 - confidence
    n = len(vals)
    
    box = np.ones((len(vals), len(vals)))
    box = box * 5
    boxlist = []
    
    for r in range(len(vals)):
        for c in range(len(vals)):
            if(r > c):
                box[r,c] = (vals[r] - vals[c]) / (r-c)
                boxlist.append((vals[r] - vals[c]) / (r-c))
    freq = 0
    #Lets caluclate frequency now
    tp = np.unique(vals, return_counts = True)
    for tpx in range(len(tp[0])):
        if(tp[1][tpx]>1):
            tp1 = tp[1][tpx]
            sev = tp1 * (tp1 - 1) * (2 * tp1 + 5)
            freq = freq + sev
        
    se = ((n * (n-1) * (2 * n + 5) - freq) / 18.0) ** 0.5
    
    no_of_vals = len(boxlist)
    
    #lets find K value
    
    k = st.norm.ppf(1-(0.05/2)) * se
    
    slope = np.median(boxlist)
    return slope, k, se

#ras = gdal.Open(inp)
#no_of_bands = ras.RasterCount               #Change this for seasonal

from glob import glob
fl = glob(inpfol + "\\*.tif")


print("\n{} files found!!!".format(len(fl)))


bigarr = []

for f in fl:
    ras = gdal.Open(f)
    arr = np.array(ras.ReadAsArray())
    
    bigarr.append(arr)
'''
for x in range(1,no_of_bands + 1):
    band = ras.GetRasterBand(x)
    bigarr.append(band.ReadAsArray())
'''    
row = ras.RasterYSize
col = ras.RasterXSize

#-Output vars------------------------------------------------------------------
z = np.zeros((row,col))
p = np.zeros((row,col))
tr_type = np.zeros((row,col))
slope = np.zeros((row,col))

print("\n")
        
#------------------------------------------------------------------------------

for r in tqdm(range(row)):
    for c in range(col):
        #sys.stdout.write("\r r = " + str(r) + " c = " + str(c) + "         ")
        #sys.stdout.flush()
        if(bigarr[0][r][c] != -999):
                vals = []
                
                for x in range(len(bigarr)):
                    vals.append(bigarr[x][r][c])
                
                
                
                z[r,c] = mann_kendall(vals,confidence)[0]
                p[r,c] = mann_kendall(vals,confidence)[1]
                tr_type[r,c] = mann_kendall(vals,confidence)[2]
                slope[r,c] = sen_slope(vals,confidence)[0]
               
            
#----------Export Outputs------------------------------------------------------        

rows = ras.RasterYSize
cols = ras.RasterXSize
geotransform = ras.GetGeoTransform()
originX = geotransform[0]
originY = geotransform[3]
pixelw = geotransform[1]
pixelh = geotransform[5]

fn = outfol + "\\Z.tif"
    
    
driver = gdal.GetDriverByName("GTiff")
outRaster = driver.Create(fn, cols, rows, 1, gdal.GDT_Float32,)
outRaster.SetGeoTransform((originX,pixelw,0,originY,0,pixelh))
outRaster.GetRasterBand(1).WriteArray(z)
outRasterSRS = osr.SpatialReference()
outRasterSRS.ImportFromWkt(ras.GetProjectionRef())
outRaster.SetProjection(outRasterSRS.ExportToWkt())
outRaster.FlushCache()
            
        
fn = outfol + "\\P.tif"
    
    
driver = gdal.GetDriverByName("GTiff")
outRaster = driver.Create(fn, cols, rows, 1, gdal.GDT_Float32,)
outRaster.SetGeoTransform((originX,pixelw,0,originY,0,pixelh))
outRaster.GetRasterBand(1).WriteArray(p)
outRasterSRS = osr.SpatialReference()
outRasterSRS.ImportFromWkt(ras.GetProjectionRef())
outRaster.SetProjection(outRasterSRS.ExportToWkt())
outRaster.FlushCache()            
            
            
fn = outfol + "\\trend.tif"
    
    
driver = gdal.GetDriverByName("GTiff")
outRaster = driver.Create(fn, cols, rows, 1, gdal.GDT_Float32,)
outRaster.SetGeoTransform((originX,pixelw,0,originY,0,pixelh))
outRaster.GetRasterBand(1).WriteArray(tr_type)
outRasterSRS = osr.SpatialReference()
outRasterSRS.ImportFromWkt(ras.GetProjectionRef())
outRaster.SetProjection(outRasterSRS.ExportToWkt())
outRaster.FlushCache()
            
        
fn = outfol + "\\slope.tif"
    
    
driver = gdal.GetDriverByName("GTiff")
outRaster = driver.Create(fn, cols, rows, 1, gdal.GDT_Float32,)
outRaster.SetGeoTransform((originX,pixelw,0,originY,0,pixelh))
outRaster.GetRasterBand(1).WriteArray(slope)
outRasterSRS = osr.SpatialReference()
outRasterSRS.ImportFromWkt(ras.GetProjectionRef())
outRaster.SetProjection(outRasterSRS.ExportToWkt())
outRaster.FlushCache()         
