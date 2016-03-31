#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import argparse,os
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
import pandas as pd

def grep(l,s):
    return [i for i in l if s in i]

def Graylevelrunlengthmatrixfeatures(ngrl):
    '''Function calculate Gray level run length matrix based features'''
    dircontent = os.listdir('.')
    selionimg = grep(dircontent,'.sim')
    fielid = []; glrlfeature=[]
    for f in selionimg:
      sname = f[:-4]
      print('PROCESSING %s' %sname)
      Img = np.genfromtxt(f,dtype=float,delimiter=',')
      if mask == 'drug':
         print('<--- Using drug mask -->')
         Mask = np.genfromtxt(sname + '_drug.msk',dtype=float,delimiter=',')
       if mask == 'mim':
         print('<--- Using MIM tissue mask -->')
         Mask = np.genfromtxt(sname + '_mim.msk',dtype=float,delimiter=',')
       if mask == 'tic':
         print('<--- Using TIC tissue mask -->')
         Mask = np.genfromtxt(sname + '_mim.msk',dtype=float,delimiter=',')
         
      ## rescaling to the desired number of gray levels
      if (ngrl != 0):
        m = ngrl/Img.max()
        scaledImg = Img*m
        binnedImg = np.rint(scaledImg)
        Img = (binnedImg + 1)  
      else:
        Img = (Img +1)   
      tissue = np.multiply(Img,Mask) 
      tissue = pd.DataFrame(tissue)
      rdf = com.convert_to_r_dataframe(tissue)
      ro.globalenv['tissue'] = rdf
      ro.r('tissue <- as.matrix(tissue)')
      ro.r('library(radiomics)')
      ro.r('glrlmatrix        <- glrlm(tissue)')
      ro.r('glrlmatrix[0,]    <- 0')                           ### Assign zero value to first row which belongs to mask region
      ro.r('glrlfeature       <- array(NA,dim=c(11,1))')
      ro.r('glrlfeature[1,1]  <- glrlm_GLN(glrlmatrix)')
      ro.r('glrlfeature[2,1]  <- glrlm_HGLRE(glrlmatrix)')
      ro.r('glrlfeature[3,1]  <- glrlm_LRE(glrlmatrix)')
      ro.r('glrlfeature[4,1]  <- glrlm_LRHGLE(glrlmatrix)')
      ro.r('glrlfeature[5,1]  <- glrlm_LRLGLE(glrlmatrix)')
      ro.r('glrlfeature[6,1]  <- glrlm_LGLRE(glrlmatrix)')
      ro.r('glrlfeature[7,1]  <- glrlm_RLN(glrlmatrix)')
      ro.r('glrlfeature[8,1]  <- glrlm_RP(glrlmatrix)')
      ro.r('glrlfeature[9,1]  <- glrlm_SRE(glrlmatrix)')
      ro.r('glrlfeature[10,1] <- glrlm_SRHGLE(glrlmatrix)')
      ro.r('glrlfeature[11,1] <- glrlm_SRLGLE(glrlmatrix)')
      glr = ro.r.matrix(ro.r('glrlfeature'))
      glr = np.array(glr)
      glrlfeature.append(glr.transpose())
      fielid.append(sname)
    glrlfeature = np.array(glrlfeature)    
    output = pd.DataFrame(glrlfeature.reshape(glrlfeature.shape[0],glrlfeature.shape[2]),columns = ["GLN","HGLRE","LRE","LRHGLE","LRLGLE","LGLRE","RLN","RP","SRE","SRHGLE","SRLGLE"])
    output.insert(0,'Id',fielid)
    output.to_csv("GLRLM_features.csv",delimiter=",")


def main():    
    parser = argparse.ArgumentParser(description="Size-zone matrix based features calculation")
    parser.add_argument('-n', dest="nugrlvl", type=int, default=32, help="Define nu of gray-level for Size Zone Matrix Calculations. (Default: 32)")
    parser.add_argument('-mask',dest = "mask",type = str, default='drug',help = "Mask used for GLRLM calculations. Available options 'drug', 'tic', 'mim'. (Default: 'drug')")
    args = parser.parse_args()
    Graylevelrunlengthmatrixfeatures(ngrl=args.nugrlvl,mask = args.mask)

        
if __name__ == '__main__':
   main()
