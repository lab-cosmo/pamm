#!/usr/bin/python
from numpy import *
from scipy.special import sph_harm
import sys
 
natoms = 38
rc = float(sys.argv[2]) 
sig = 2
lmax = 6 

ixyz = open(sys.argv[1])
 
x = zeros([natoms,3])

frame = 0
while True:
  line = ixyz.readline()
  if line=='': break
  natoms = int(line)
  
  # skip comment
  comment = ixyz.readline()
  
  for i in xrange(natoms):
    # loop over positions
    xyz = ixyz.readline()
    xyz = xyz.split()
    x[i,0] = float(xyz[1])
    x[i,1] = float(xyz[2])
    x[i,2] = float(xyz[3])

  # compute fermi function and
  # modified steinhardt parameter
  for i in xrange(natoms):
    rij = x[i]-x[i != arange(natoms)]
    r = linalg.norm(rij,axis=1)
    phi = arccos(rij[:,2]/r)
    theta = arctan2(rij[:,1],rij[:,0]) + pi
    
    # fermi function
    f = 1./(exp((r-rc)/sig)+1.)
    # sum of all fermi contributions
    sf = sum(f)
  
    # steinhardt order parameter
    q = zeros(lmax)
    for l in xrange(4,lmax+1):
      for m in xrange(-l,l+1):
        qlm = 1./sf * sum(sph_harm(m,l,theta,phi) * f)
        q[l-1] = q[l-1] + abs(qlm)**2
      q[l-1] = sqrt(q[l-1] * 4.*pi/(2.*l + 1.))

    q[1]=sf
    print ' '.join(map(str,q[1:6:2]))
  
#  print 'frame %d' % frame
  frame = frame + 1

ixyz.close()


