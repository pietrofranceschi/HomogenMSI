#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from scipy import misc
from scipy import ndimage
from PIL import Image
import matplotlib.pyplot as plt
import os, struct
import numpy as np
from skimage import filter
import argparse
import itertools
### Functions to read the Analyze 7.5 format
def readAnalyzeheader(filename):
   '''Function reading Analyze header file
   
   The header contains informations regardin gthe number of 
   points and other dataset characteristics.

   Args:
   filename : The name of the hdr file

   Value:
   nx : number of pixels in th x direction
   ny : number of pixel in the y direction
   xd : step in x 
   xy : step in y
   '''
   size = os.stat(filename).st_size
   if((size != 348) & (size != 384)):
       sys.exit("%s file is not Anlayze format.", filename)
   dimensions = np.zeros(8)
   pixdim = np.zeros(8)
   f = open(filename,'rb')
   f.seek(38)
   regular = f.read(1)
   if(regular != 'r'):
       sys.exit("wrong file format")
   f.seek(40)
   for i in range(8):
       dimensions[i] = struct.unpack('H',f.read(2))[0]
   nx = dimensions[2]  # x size of the image (number of columns)
   ny = dimensions[3]  # y size of the image (number of rows)
   f.seek(70)
   datatype = struct.unpack('H',f.read(2))[0]
   bitpix = struct.unpack('H',f.read(2))[0]
   if(datatype == 2 or datatype == 4 or datatype == 8 ):
       what = 'integer'
   elif (datatype == 16 or datatype == 32 or datatype == 64 ):
       what = 'double'
   else:
       what = 'raw'
   signed = (datatype == '2') 
   size = (bitpix/8)
   f.seek(76) 
   for i in range(8):
       pixdim[i] = struct.unpack('f',f.read(4))[0]
   xd = pixdim[1]
   yd = pixdim[2]
   return(int(nx),int(ny), xd, yd, signed, size, what) 

## read t2m

def readAnalyzet2m(filename):
    '''Function reading t2m mass file '''
    totalBytes = os.stat(filename).st_size
    bytes = 4
    mass = np.zeros(totalBytes/bytes)
    endian = 'f'
    with open(filename,'rb') as f:
        for i in xrange(len(mass)):
            mass[i] = struct.unpack(endian,f.read(bytes))[0]
    return np.array(mass)  

## read the image

def readAnalyzeImage(filename,t2m,hdr,massrange = []):
    '''
    This function  interacts with readAnalyzet2m and readAnalyzeheader
    and extract 2d map showing the integrated intensity in the data
    '''
    if (massrange == []):
        mass1, mass2 = min(t2m), max(t2m)
    else:
        mass1, mass2 = massrange[0], massrange[1]
    id = (t2m >= mass1) & (t2m <= mass2)
    bytes, endian = 2, 'H'
    datamatrix = np.zeros([hdr[0],hdr[1]]) ## the matrix to store the image
    with open(filename,'rb') as f:
        for c,r in itertools.product(range(hdr[1]),range(hdr[0])):
            spcarr = []
            for i in range(t2m.size):
                if not id[i]:
                    f.seek(bytes,1)
                    continue 
                test = struct.unpack(endian,f.read(bytes))[0]
                spcarr.append(test)
            datamatrix[r,c] = sum(spcarr)
    return(datamatrix)

def getidmaxIntensity(filename,t2m):
    spectra = []
    f = open(filename,'rb')
    f.seek(t2m.size * int(20) * 2)
    for i in xrange(t2m.size):
        test = struct.unpack('H',f.read(2))[0]
        spectra.append(test)
    spectra = np.array(spectra)    
    nu =  np.where(spectra == max(spectra))[0][0]
    return(nu)

def grep(l,s):
    return [i for i in l if s in i]

## batch process stuff ...
def processMSIBatch(massrange=[], mfiltrad = 3, msk_tic = False, msk_mim = False,massrange_tissue = False, massrange_std = False):
    '''
    Process in batch mode all the analyze files contained in a FOlder
    '''
   
    dircontent = os.listdir('.')  ## dirrctory 
    hdrs = grep(dircontent,'.hdr')
    for f in hdrs:
        myheader = readAnalyzeheader(f)
        sname = f[:-4]      ## get rid of the extension
        print('PROCESSING %s : -mzmin %.3f -mzmax %.3f' % (sname,massrange[0],massrange[1]))
        
        mass = readAnalyzet2m(sname + '.t2m')
        ## extract and save the drug ion image
        if massrange_tissue:
               MSImatrix_drug = readAnalyzeImage(sname + '.img',mass,myheader,massrange[0:2]) 
               MSImatrix_tissue = readAnalyzeImage(sname + '.img',mass,myheader,massrange[2:4])              
               MSImatrix_drug = np.sqrt(MSImatrix_drug); MSImatrix_tissue = np.sqrt(MSImatrix_tissue); 
               MSImatrix_drug = MSImatrix_drug +1   # To make necrosis region different from background zero value

############## Dividing drug and tissue image with standard              
               if massrange_std:
		   MSImatrix_std = readAnalyzeImage(sname + '.img',mass,myheader,massrange_std) 
		   MSImatrix_std = np.sqrt(MSImatrix_std)
		   MSImatrix_std = MSImatrix_std + 1
                   MSImatrix_drug = MSImatrix_drug / MSImatrix_std
                   MSImatrix_tissue = MSImatrix_tissue / MSImatrix_std
                             
############## Object detection based on tissue mass  
               print('<--- Extracting and saving tissue object -->')        
               MSImatrix_tissue = ndimage.median_filter(MSImatrix_tissue,mfiltrad)
               nu = []
               nr = MSImatrix_tissue.shape[0]-1
               nc = MSImatrix_tissue.shape[1]-1
               nu.append(MSImatrix_tissue[0,0:4])
               nu.append(MSImatrix_tissue[0,(nc-5):(nc-1)])
               nu.append(MSImatrix_tissue[nr,0:4])
               nu.append(MSImatrix_tissue[nr,(nc-5):(nc-1)])
               flattened_list = [y for x in nu for y in x]
               mask_tissue = MSImatrix_tissue > np.median(flattened_list)
               mask_tissue = mask_tissue.astype(int)
               
############## Create drug mask   
               print('<--- Exracting and Saving Drug Image -->')
               mask_drug = MSImatrix_drug * mask_tissue
               mask_drug = ndimage.median_filter(mask_drug,mfiltrad)
               val = filter.threshold_otsu(mask_drug)
               mask_tissue = mask_drug> val
                              
               np.savetxt(sname + '.sim',mask_drug,fmt='%-3.2f',delimiter=',')
               np.savetxt(sname + '_drug.msk',mask_tissue,fmt ='%-3.2f',delimiter=',')
               fig = plt.figure(1)
               ax1 = plt.subplot(211)
               plt.title(sname)
               plt.imshow(mask_drug,interpolation='None')
               plt.colorbar()
               ax2 = plt.subplot(212)
               plt.title('drug mask')
               plt.imshow(mask_tissue,interpolation='None',cmap='Greys')
               plt.colorbar()
               plt.savefig(sname + '.jpg',bbox_inches='tight')
               plt.close()

        else:
               print('<--- Exracting and Saving Drug Image -->')
               MSImatrix = readAnalyzeImage(sname + '.img',mass,myheader,massrange) 
               
############## Dividing drug and tissue image with standard              
               if massrange_std:
		   MSImatrix_std = readAnalyzeImage(sname + '.img',mass,myheader,massrange_std) 
		   MSImatrix_std = np.sqrt(MSImatrix_std)
		   MSImatrix_std = MSImatrix_std + 1
                   MSImatrix = MSImatrix / MSImatrix_std
                   
               MSImatrix = ndimage.median_filter(MSImatrix,mfiltrad)

############## Create drug mask
               print('<--- Exracting and Saving Drug Mask -->')
               val = filter.threshold_otsu(np.sqrt(MSImatrix))
               mask_drug = np.sqrt(MSImatrix) > val
               mask_drug = mask_drug.astype(int)        

               np.savetxt(sname + '.sim',MSImatrix,fmt='%-3.2f',delimiter=',')
               np.savetxt(sname + '_drug.msk',mask_drug,fmt ='%-3.2f',delimiter=',')
               fig = plt.figure(1)
               ax1 = plt.subplot(211)
               plt.title(sname)
               plt.imshow(MSImatrix,interpolation='None')
               plt.colorbar()
               ax2 = plt.subplot(212)
               plt.title('drug mask')
               plt.imshow(mask_drug,interpolation='None',cmap='Greys')
               plt.colorbar()
               plt.savefig(sname + '.jpg',bbox_inches='tight')
               plt.close()


####### Create TIC image mask
        if msk_tic:
        	print('<--- Creating and Saving TIC based Tissue Mask -->')
        	MSImatrix_tic = readAnalyzeImage(sname + '.img',mass,myheader,massrange = []) 
        	MSImatrix_tic = ndimage.median_filter(MSImatrix_tic,mfiltrad)
        	val = filters.threshold_otsu(np.sqrt(MSImatrix_tic))
        	mask_tic = np.sqrt(MSImatrix_tic) > val
        	mask_tic = mask_tic.astype(int)        
        	np.savetxt(sname + '_tic.msk',mask_tic,fmt='%-3.2f',delimiter=',')
        	## Save the mask
        	plt.title(sname + ' tissue_tic')
        	plt.imshow(mask_tic,interpolation='None',cmap='Greys')
        	plt.savefig(sname + '_ticmsk'+'.jpg',bbox_inches='tight')
        	plt.close()

######## Create background mask with the maximum intensity mass	
        if msk_mim:
        	print('<--- Creating and Saving MIM based Tissue Mask -->')
        	ind = getidmaxIntensity(sname + '.img',mass)
        	mass1 = mass[ind]-0.8; mass2 = mass[ind]+0.8
        	MSImatrix_mim = readAnalyzeImage(sname + '.img',mass,myheader,[mass1,mass2])
        	## otsu filtering
        	MSImatrix_mim = ndimage.median_filter(MSImatrix_mim,mfiltrad)
        	val = filters.threshold_otsu(np.sqrt(MSImatrix_mim))
        	mask_mim = np.sqrt(MSImatrix_mim) > val
        	mask_mim = mask_mim.astype(int)
        	mask_mim = np.absolute(mask_mim-1)
        	np.savetxt(sname + '_mim.msk',mask_mim,fmt='%-3.2f',delimiter=',')
        	plt.title(sname + 'tissue_mim')
        	plt.imshow(mask_mim,interpolation='None',cmap='Greys')
        	plt.savefig(sname + '_mimmsk'+'.jpg',bbox_inches='tight')
        	plt.close()
        

def main():    
    parser = argparse.ArgumentParser(description="Create an extracted ion image and its corresponding binary mask")
    parser.add_argument('-mz',dest = "massrange",type = float, nargs ='+', default=[],help = "Desired m/z range ")
    parser.add_argument('-mz_tissue',dest = "massrange_tissue",type = float, nargs ='+', default= False,help = "m/z range correspond to tissue")
    parser.add_argument('-mz_std',dest = "massrange_std",type = float, nargs ='+', default= False,help = "m/z range for standard")
    parser.add_argument('-fr',dest = "mfiltrad",type = int, default=3,help = "Radius of the median filter")
    parser.add_argument('-tic',dest = "msk_tic",action='store_true',help = "TIC based tissue identification", default = False)
    parser.add_argument('-mim',dest = "msk_mim",action='store_true',help = "Maximum Intensity Ion tissue identification", default = False)
    args = parser.parse_args()
    processMSIBatch(massrange = args.massrange,massrange_std = args.massrange_std ,massrange_tissue = args.massrange_tissue ,mfiltrad = args.mfiltrad, msk_tic = args.msk_tic, msk_mim=args.msk_mim)
                    
if __name__ == '__main__':
   main()








