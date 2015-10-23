import numpy as np
import argparse
import cv2,os
import pandas as pd

def grep(l,s):
    return [i for i in l if s in i]

def Shapefactorfeatures():
    '''Function calculate shape factors statistics for tissue mask images'''
    dircontent = os.listdir('.')
    selmaskimg = grep(dircontent,'_drug.msk')
    stat = []
    for f in selmaskimg:
        sname = f[:-9];areat=0;perimetert = 0
        print('PROCESSING %s' %sname)
        mask_drug = np.genfromtxt(f,dtype=float,delimiter=',')
        shapefactors = {}
        # Analysis using bkg mask image to get all global parameters
        contours,_ = cv2.findContours(np.uint8(mask_drug),cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)   #Finding contours
        nregions = len(contours)  ## number of non conected regions
        areas = [cv2.contourArea(cnt) for cnt in contours]
        area = sum(areas)
        perimeters = [cv2.arcLength(cnt,True) for cnt in contours]
        perimeter = sum(perimeters)
        shapefactors['Id'] = sname
        shapefactors['nregions'] = nregions
        shapefactors['global_area'] = area
        shapefactors['global_perimeter'] = perimeter
        shapefactors['global_circularity'] = ((4*np.pi*area)/(perimeter**2))
        stat.append(shapefactors)
    output = pd.DataFrame(stat)
    output.to_csv("SB_features.csv",sep=',')
        
def main():    
    parser = argparse.ArgumentParser(description="Calcuates the shape factors for all drug mask images")
    args = parser.parse_args()
    Shapefactorfeatures()

        
if __name__ == '__main__':
   main()