import numpy as np
import os
from scipy import stats
import argparse
import pandas as pd

def grep(l,s):
    return [i for i in l if s in i]


def Firstorderstatistics(mask,ngl):
    '''Function calculate first order statistics the tissue image'''
    dircontent = os.listdir('.')  ## dirrctory content
    selionimg = grep(dircontent,'.sim')
    stat = []
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
            Mask = np.genfromtxt(sname + '_tic.msk',dtype=float,delimiter=',')
        ## rescaling to the desired number of gray levels
        m = ngrl/Img.max()
        scaledImg = Img*m
        binnedImg = np.rint(scaledImg)
        Img = (binnedImg + 1)  

        tissueImg = np.multiply(Img,Mask)
        tissue = tissueImg[np.where(tissueImg>0)]
        FOS = {}
        FOS['Id'] = sname
        FOS['sum'] = np.sum(tissue)            
        FOS['pixel_counts'] = np.sum(tissue > 0.001) 
        FOS['mean_intensity'] = np.mean(tissue)
        FOS['median_intensity'] = np.median(tissue)
        FOS['standard_deviation'] = np.std(tissue)
        FOS['variance'] = np.var(tissue)
        FOS['kurtosis'] = stats.kurtosis(tissue.flatten())
        FOS['skewness'] = stats.skew(tissue.flatten())
        stat.append(FOS)
    output = pd.DataFrame(stat)
    output.to_csv("FOS_features.csv", sep=',')



def main():    
    parser = argparse.ArgumentParser(description="Calculates the first order statistics of a series of .sim images. The outputs are stored as a csv")
    parser.add_argument('-mask',dest = "mask",type = str, default='drug',help = "Mask used for FOS calculations. Available options 'drug', 'tic', 'mim'. (Default: 'drug')")
    parser.add_argument('-n', dest="nugrlvl", type=int, default=32, help="Define number of gray-level for FOS calculation (Default: 32)")
    args = parser.parse_args()
    Firstorderstatistics(mask = args.mask,ngrl = args.nugrlvl)
        
if __name__ == '__main__':
   main()
