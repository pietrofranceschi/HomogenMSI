#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cmath,os
import pandas as pd
import numpy as np
import argparse
from skimage.feature import greycomatrix, greycoprops

def _entropy(p):
    ''' Function calcuate entropy feature'''
    p = p.ravel()
    return -np.dot(np.log2(p+(p==0)),p)

def gray_features(p):
    '''
    Function calculate all co-occurence matrix based features
    Input p: gray level co-occurence matrix

    Output: List of the 13 GCLM features

    '''
    feats = np.zeros(13,np.double)
    maxv = len(p)
    k = np.arange(maxv)
    k2 = k**2
    tk = np.arange(2*maxv)
    tk2 = tk**2
    i,j = np.mgrid[:maxv,:maxv]
    ij = i*j
    i_j2_p1 = (i - j)**2
    i_j2_p1 += 1
    i_j2_p1 = 1. / i_j2_p1
    i_j2_p1 = i_j2_p1.ravel()
    px_plus_y = np.empty(2*maxv, np.double)
    px_minus_y = np.empty(maxv, np.double)
    pravel = p.ravel()
    px = p.sum(0)
    py = p.sum(1)
    ux = np.dot(px, k)
    uy = np.dot(py, k)
    vx = np.dot(px, k2) - ux**2
    vy = np.dot(py, k2) - uy**2
    sx = np.sqrt(vx)
    sy = np.sqrt(vy)
    px_plus_y = np.zeros(shape=(2*p.shape[0] ))
    px_minus_y = np.zeros(shape=(p.shape[0]))
    for i in range(p.shape[0]):
       for j in range(p.shape[0]):
           p_ij = p[i,j]
           px_plus_y[i+j] += p_ij
           px_minus_y[np.abs(i-j)] += p_ij
    feats[0] = np.sqrt(np.dot(pravel, pravel))                        # Energy
    feats[1] = np.dot(k2, px_minus_y)                                 # Contrast
    if sx == 0. or sy == 0.:
       feats[2] = 1.
    else:
       feats[2] = (1. / sx / sy) * (np.dot(ij.ravel(), pravel) - ux * uy) # Correlation
    feats[3] = vx                                                     #Sum of Squares: Variance     
    feats[4] = np.dot(i_j2_p1, pravel)                                # Inverse of Difference Moment
    feats[5] = np.dot(tk, px_plus_y)                                  # Sum Average
    feats[7] = _entropy(px_plus_y)                                    # Sum Entropy
    feats[6] = ((tk-feats[7])**2*px_plus_y).sum()                     # Sum Variance
    feats[8] = _entropy(pravel)                                       # Entropy
    feats[9] = px_minus_y.var()                                       # Difference Variance
    feats[10] = _entropy(px_minus_y)                                  # Difference Entropy
    HX = _entropy(px)
    HY = _entropy(py)
    crosspxpy = np.outer(px,py)
    crosspxpy += (crosspxpy == 0) 
    crosspxpy = crosspxpy.ravel()
    HXY1 = -np.dot(pravel, np.log2(crosspxpy))
    HXY2 = _entropy(crosspxpy)
    if max(HX, HY) == 0.:
       feats[11] = (feats[8]-HXY1)                                    # Information Measure of Correlation 1                    
    else:
       feats[11] = (feats[8]-HXY1)/max(HX,HY)
    feats[12] = np.sqrt(max(0,1 - np.exp( -2. * (HXY2 - feats[8]))))  # Information Measure of Correlation 2
    return feats


def grep(l,s):
    return [i for i in l if s in i]

def GLCMfeatures(d,ngrl,mask):
    '''
    Function for feature calculation in all possible direction, for given distance parameter
    '''
    
    dircontent = os.listdir('.')
    selionimg = grep(dircontent,'.sim')
    stat = [];
    #default angle parameter value for co-occurence matrix calculation
    a = [0, np.pi/4, np.pi/2, 3*np.pi/4]     
    
    #Start calculation for individual dataset'''
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
        ## calculate the masked tissue
        tissue = np.multiply(Img,Mask)
        GLCM = {}
        #GLCM calculation for individual distance value in all four directions
        feat_dist=[]         
        for i in d:                                     
            feat_angle = []
            for j in a:                               # for loop on angle parameter
                g = greycomatrix(tissue,[i],[j])
                g = g.reshape(g.shape[0],g.shape[1])
                g[0,0] = 0                                         # providing zero value to gray level corresponding to background region
                g = g.astype(np.float)
                g = (g/np.sum(g))                                  # Normalization of co-occurrence matrix
                feat_angle.append(gray_features(g))                # Extract feature values for co-coccurrence matrix g
            feat_angle = np.array(feat_angle) 
            feat_dist.append(np.mean(feat_angle,axis=0))            # Averaging over all four directions  
            ## collect the GCLM features for all the images
        GLCM['Id']                      = sname
        GLCM['energy']                  = (np.mean(feat_dist,axis=0))[0]
        GLCM['contrast']                = (np.mean(feat_dist,axis=0))[1]
        GLCM['correlation']             = (np.mean(feat_dist,axis=0))[2]
        GLCM['sum_squares_Variance']    = (np.mean(feat_dist,axis=0))[3]
        GLCM['homogeneity']             = (np.mean(feat_dist,axis=0))[4]
        GLCM['sum_average']             = (np.mean(feat_dist,axis=0))[5]
        GLCM['sum_variance']            = (np.mean(feat_dist,axis=0))[6]
        GLCM['sum_entropy']             = (np.mean(feat_dist,axis=0))[7]
        GLCM['entropy']                 = (np.mean(feat_dist,axis=0))[8]
        GLCM['diff_variance']           = (np.mean(feat_dist,axis=0))[9]
        GLCM['diff_entropy']            = (np.mean(feat_dist,axis=0))[10]
        GLCM['iMC1']                    = (np.mean(feat_dist,axis=0))[11]
        GLCM['iMC2']                    = (np.mean(feat_dist,axis=0))[12]
        stat.append(GLCM)
    output = pd.DataFrame(stat)
    output.to_csv("GLCM_feaures.csv", sep=',')


def main():    
    parser = argparse.ArgumentParser(description="Calculates the first order statistics of a series of .sim images. The outputs are stored as a csv")
    parser.add_argument('-d',dest="distance",type = int,nargs='+', default=[1], help="Distance parameter for co-occurrence matrix calculation")
    parser.add_argument('-n', dest="nugrlvl", type=int, default=32, help="Define number of gray-level for co-occurrence matrix calculation (Default: 32)")
    parser.add_argument('-mask',dest = "mask",type = str, default='drug',help = "Mask used for GCLM calculations. Available options 'drug', 'tic', 'mim'. (Default: 'drug')")
    args = parser.parse_args()
    GLCMfeatures(d = args.distance,ngrl = args.nugrlvl, mask = args.mask)
    #print(args)
           
if __name__ == '__main__':
   main()