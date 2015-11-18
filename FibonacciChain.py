'''
Written by Konrad on
2015-05-27 Wed 05:56 AM
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
import matplotlib.cm as cm
import csv
from scipy.optimize import curve_fit
#print help(curve_fit)

path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not path in sys.path:
    sys.path.insert(1, path)
    del path

import functionsZH as ZH
reload(ZH)

import phys_const as pc
import K39_const as K39

tau = 0.5*(1.+np.sqrt(5.))

N = 300

start = 0
stop = 20


a = np.arange(start,stop+1,1).astype(np.int)
print a

def integerPart(x):
    return np.floor((x+1.)/tau).astype(np.int)

print integerPart(a)


def fibonacciDistance(n):
    '''Returns distance from origin to an atom at the n-th position in a fibonacci chain with segment lengths S = 1 and L = tau. First term: add one short segment in each step, second term: add an extra (L-S) if the segment is long.'''
    return n + 1./tau*np.floor((n+1)/tau)


def kBragg(p,q):
    return 2.*np.pi*tau**2/(tau**2+1.)*(p/tau+q)

def structureFactor(p,q):
    return np.sin(np.pi*tau/(tau**2+1.)*(tau*p-q))/(np.pi*tau/(1.+tau**2)*(tau*p - q))*np.exp(1j*np.pi*(tau - 2.)/(tau + 2.)*(tau*p-q))

def getAllCombinations(a):
    '''creates a (len(a),2)-dim array with all possible combinations of a with itself. Takes only 1D arrays.'''
    out = np.zeros((len(a)**2,2))
    out[:,0] = np.repeat(a, len(a))
    for i0 in np.arange(0,len(a)):
        out[len(a)*i0:len(a)*(i0+1), 1] = a[:]
    return out



momentumTransferArray = np.array([])
structureFactorArray = np.array([])

#for p in np.arange(0,N,1):
#    for q in np.arange(0, N,1):
#        momentumTransferArray = np.append(momentumTransferArray,kBragg(p,q))
#        structureFactorArray = np.append(structureFactorArray, structureFactor(p,q))


my_range = np.arange(0,N)
my_combns = getAllCombinations(my_range)

momentumTransferArray = kBragg(my_combns[:,0],my_combns[:,1])
structureFactorArray = structureFactor(my_combns[:,0],my_combns[:,1])
print structureFactorArray.shape


intensityArray =  np.abs(structureFactorArray)
#print np.amax(intensityArray)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(momentumTransferArray, intensityArray, 'b.')
plt.show()




