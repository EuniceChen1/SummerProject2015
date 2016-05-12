from __future__ import division, print_function
import numpy as np
import os
import glob
import healpy as hp
from astropy.table import Table
import argparse
import ast
import pandas as pd
#--------------------------------- LOGIC

''' 
Have multiple unprocessed galaxy catalogs. 
Use a for loop and input those catalogs to cut them one by one.
Get processed catalogs from the for loop.
Input processed catalogs into bin_data function and bin the processed catalogs.
Use a for loop in the bin_data function to bin the catalogs one by one.
Binned catalogs inputted into the main function for loop to create a full count map.

'''
#-----------------------------------

def apply_cuts(data):

    #sdss = Table.read(fnames, format='csv')

    # Part 1: Filtering the catalog with only galaxy type
    separation = data['type']   
    galaxy = (separation == 6)
    data = data[galaxy]    

    # Part 2: Choosing on clean photometry (sdss['clean'] = 1)

    cleaning = data['clean']
    clearimage = (cleaning == 1)
    data = data[clearimage]

    # Part 3: Choosing image mask

    mask = data['insideMask']
    maskmap = (mask == 0)
    data = data[maskmap]

    newdata = data['ra','dec','type','clean','insideMask']
	
    return newdata
	
def bin_data(data,field,binlist):

    bincolumn = data[field]
    bindata = ((bincolumn >= inf) & (bincolumn < sup) for (inf, sup) in binlist)
    datalist = [data[mask] for mask in bindata]

    return datalist

def main(infile, nside, RA, DEC, zname, zbins, cuts):
	
    filelist = glob.glob(infile)
    filelist.sort()
    folder = os.getcwd()
	
    # Create empty Healpix map
    full_countmap = np.zeros(hp.nside2npix(nside))
    # Get number of pixels from the map
    npix = hp.nside2npix(nside)
	
    for filename in filelist:
	# 1) Open catalogs and perform checks
	filepath = os.path.join(folder, filename)
	print (filepath) 
	# Open the catalog as a table
	header_row=["objID","ra","dec","type","clean","insideMask","z","zErr","photoErrorClass","nnCount","chisq","rnorm","bestFitTemplateID","absMagU","absMagG","absMagR","absMagI","absMagZ"]
	data = Table(np.array(pd.read_csv(filepath,skiprows=1)),header_row)
	#Check if the file is empty
	if (len(data) == 0):
	    print ("file is empty")
	    continue
	else:
	    print ("file is not empty")
	   	#Renaming columns
	    selectcols = ["col0","col1","col2","col3","col4","col5","col6"]
	    newcolnames = ["objID","ra","dec","type","clean","insideMask","z"]
	    for oldname, newname in zip(selectcols, newcolnames):
			data.rename_column(oldname,newname)
	    colnames = data.colnames
	    # Make sure that RA and DEC columns are in the catalog
	    assert (RA in colnames) and (DEC in colnames), ("Both ra and dec must" +
							"be in the catalog")

	   # Make sure that Nside is a power of 2
	    assert hp.isnsideok(nside), "nside must be a power of 2"

		
	    # 2) Apply general cuts and bin in redshift
	    # Apply general cuts to catalog
	    if cuts is not None:
		data = apply_cuts(data)
		    
	    # Apply redshift bin
	    if (zname is not None) and (zbins is not None):
		datalist = bin_data(data, zname, zbins)
		zsuffix = ["_z%.2f-%.2f" % (inf, sup) for (inf, sup) in zbins]
	    else:
		print ("For redshift binning, both zname and zbins must be given.")
		print ("Proceeding without binning.")
		datalist = [data]
		zsuffix = [""]

	    # 3) Create count maps for each redshift bin
	    # Translate radec coordinates to healpix format
	    theta = np.deg2rad(90.0 - data['dec'])
	    phi = np.deg2rad(data['ra'])
		
	    # Affect each galaxy to a healpix pixel
	    try:
		gal_hppix = hp.ang2pix(nside, theta=theta, phi=phi, nest=False)
	    except ValueError:
		print("BEWARE! Problem with RA DEC range, creating fake random map.")
		theta = np.random.uniform(0, np.pi, size=len(data))
		phi = np.random.uniform(0, 2*np.pi, size=len(data))
		gal_hppix = hp.ang2pix(nside, theta=theta, phi=phi, nest=False)
			
	    # Count number of galaxies in each pixel	
	    countmap = np.bincount(gal_hppix, minlength=npix)
	    # Make sure size of count map is the same as number of pixels
	    assert len(countmap) == npix, ("Size of count map must be the same" + 
					"as minimum length")
	    # Add counts to the originally empty map
	    full_countmap += countmap

    # Save final map
    savemap = hp.write_map("/share/splinter/visit10/sdss_full_countmap.fits", full_countmap)
    return None


if __name__ == "__main__":
	
    infile = "/share/splinter/moraes/photoz_cats/photoz_cat_*.csv"
    nside = 128
    RA = 'ra'
    DEC = 'dec'
    zname = 'z'
    zbins = [(0.,0.10),(0.11,0.20),(0.21,0.30),(0.31,0.40),(0.41,0.50),(0.51,0.60),(0.61,0.70),(0.71,0.80),(0.81,0.90),(0.91,1.00)]
    cuts = None
       
    main(infile, nside, RA, DEC, zname, zbins, cuts)
