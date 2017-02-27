# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 13:57:19 2016

@author: Rebecca
"""
from PIL import Image
from pylab import *
import numpy as np
from scipy.ndimage import filters
import cv2
import matplotlib.pyplot as plt

I = np.array(Image.open('star.jpg').convert('L'))
gray() 

wZ = 0.43
wG = 0.14
wD = 0.43

fG = np.zeros(I.shape)
fZ = np.zeros(I.shape)

Il = cv2.Laplacian(I,cv2.CV_64F)
Il = (Il != 0).astype(int)
imshow(Image.fromarray(Il))
 
x = np.zeros(I.shape)
filters.sobel(I,1,x)
y = np.zeros(I.shape)
filters.sobel(I,0,y)
G = np.sqrt(x**2+y**2)
fG = 1-(G/G.max()) 

def fDfunc(p,q):
    px = p[0]
    py = p[1]
    
    qx = q[0]
    qy = q[1]
    
    Dp = (norm(y[px][py]),norm(-x[px][py]))
    Dq = (norm(y[qx][qy]),norm(-x[qx][qy]))
    
    #L(p, q) = { q – p; if D(p) · (q – p) ≥ 0,
          #      p – q; if D(p) · (q – p) < 0
    L = abs(np.subtract(p,q))
    dp = dot(Dp,L)
    dq = dot(Dq,L)
    fD = (1/pi)*( (1/cos(dp)) + (1/cos(dq)) )
    #print(fD)
    return fD
    
def cost(p,q):
    #l(p,q) = ωZ · fZ(q) + ωD · fD (p,q) + ωG · fG(q)
    l = (wZ*fZ[q[0]][q[1]]) + (wD*fDfunc(p,q)) + (wG*fG[q[0]][q[1]])
    return l

def dijkstras(s,e):
    #Initialization of data structures
    sx = s[0]
    sy = s[1]    
    
    L = [(0,0)]
    totCost = [100.0]
    p = []
    expanded = np.zeros(I.shape, dtype=bool)
    totCost.append(0.0)
    L.append(s)
    p.append(s)
    q = s
         
    
    
    while(L and (not all(q==e))):
    
        q = L[totCost.index(min(totCost))]
    
        qx = q[0]
        qy = q[1]        
        
        expanded[qx][qy] = True
        
        N = [(qx-1,qy-1),(qx,qy-1),(qx+1,qy-1),
         (qx-1,qy),            (qx+1,qy),
         (qx-1,qy+1),(qx,qy+1),(qx+1,qy+1) ]        
        
        for r in N:
            #comment out to stop at error
            #if(r[0]<I.shape[0] and r[1]< I.shape[1]): 
                if(not expanded[r[0]][r[1]]):
                    temp = totCost[L.index(q)]+cost(q,r) 
                    if (r in L and temp < totCost[L.index(r)]):
                        ind = L.index(r)
                        L.remove(r)
                        totCost.remove(totCost[ind])
                        p.remove(r)
                    else:
                        L.append(r)
                        totCost.append(cost(q,r)+temp)
                        p.append(r)
                        
                        #uncomment to show path taken
                        #plt.plot(r[1],r[0], 'r.')
                        #plt.show()
                                          
        L.remove(q)
        totCost.remove(min(totCost))
    return p
    
imshow(I)
inpu = ginput(2)    
p = dijkstras((inpu[0][1],inpu[0][0]),(inpu[1][1],inpu[1][0]))
print(p)
x_val = [x[0] for x in p]
y_val = [x[1] for x in p]
plt.plot(x_val,y_val, 'r-')
plt.show()
