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

def Sizezonematrixfeatures(ngrl,mask):
    '''Function calculate size zone matrix based features'''
    dircontent = os.listdir('.')
    selionimg = grep(dircontent,'.sim')
    fielid = []; szmfeature=[]
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
        Img = np.sqrt(Img)
        Img = np.rint(Img)
        Img = (Img +1)   
      tissue = np.multiply(Img,Mask) 
      tissue = pd.DataFrame(tissue)
      rdf = com.convert_to_r_dataframe(tissue)
      ro.globalenv['tissue'] = rdf
      ro.r('tissue <- as.matrix(tissue)')
      ro.r('library(radiomics)')
      ro.r('szmatrix <- glszm(tissue)')
      ro.r('szmatrix[0,] <- 0')                           ### Assign zero value to first row which belongs to mask region
      ro.r('szmfeature <- array(NA,dim=c(11,1))')
      ro.r('szmfeature[1,1] <- glszm_SAE(szmatrix)')
      ro.r('szmfeature[2,1] <- glszm_LAE(szmatrix)')
      ro.r('szmfeature[3,1] <- glszm_IV(szmatrix)')
      ro.r('szmfeature[4,1] <- glszm_HILAE(szmatrix)')
      ro.r('szmfeature[5,1] <- glszm_LILAE(szmatrix)')
      ro.r('szmfeature[6,1] <- glszm_HISAE(szmatrix)')
      ro.r('szmfeature[7,1] <- glszm_LISAE(szmatrix)')
      ro.r('szmfeature[8,1] <- glszm_HIE(szmatrix)')
      ro.r('szmfeature[9,1] <- glszm_LIE(szmatrix)')
      ro.r('szmfeature[10,1] <- glszm_ZP(szmatrix)')
      ro.r('szmfeature[11,1] <- glszm_SZV(szmatrix)')
      szm = ro.r.matrix(ro.r('szmfeature'))
      szm = np.array(szm)
      szmfeature.append(szm.transpose())
      fielid.append(sname)
    szmfeature = np.array(szmfeature)    
    output = pd.DataFrame(szmfeature.reshape(szmfeature.shape[0],szmfeature.shape[2]),columns = ["sae","lae","iv","szv","zp","lie","hie","lisae","hisae","lilae","hilae"])
    output['Id'] = fielid
    output.to_csv("SZM_features.csv",delimiter=",")


def main():    
    parser = argparse.ArgumentParser(description="Size-zone matrix based features calculation")
    parser.add_argument('-n', dest="nugrlvl", type=int, default=32, help="Define nu of gray-level for Size Zone Matrix Calculations. (Default: 32)")
    parser.add_argument('-mask',dest = "mask",type = str, default='drug',help = "Mask used for SZM calculations. Available options 'drug', 'tic', 'mim'. (Default: 'drug')")
    args = parser.parse_args()
    Sizezonematrixfeatures(ngrl=args.nugrlvl,mask = args.mask)

        
if __name__ == '__main__':
   main()