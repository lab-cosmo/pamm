# multivariate bayesian bandwidth estimation
from numpy import *
import pylab as P
from matplotlib.colors import LogNorm
import os

def mm(x,n=-1,verbose=False):
   "creates grid points using minmax criterion"
   # if n not set
   if n < 0: n = int(len(x)**(1.0/2.0))
   x = array(x)
   d = 1
   if (x.ndim > 1): d = x.shape[1]
   y = empty((n,d))
   x.reshape(-1,d)
   
   dmin = array((len(x))*[float('inf')])
   
   # choosing a random point
   y[0,:] = x[random.randint(0,size(x,0)),:]
   for i in range(1,n):
     if verbose and ((i+1)%10 is 0): 
       print "minmax:  %s/%s" % (i+1,n)
     dx = linalg.norm(x-y[i-1,:],axis=1)
     dmin[dmin > dx] = dx[dmin > dx]
     y[i,:] = x[argmax(dmin),:]
   return y

def oracle(Q,N):
  "oracle shrinkage of matrix"
  
  if len(Q.shape) < 2: D = 1.
  else: D = float(Q.shape[1])
  
  # apply oracle approximating shrinkage alogorithm on local Q
  a = ( (1.-2./D) * trace(Q**2) \
         + trace(Q)**2 ) / ( (N + 1. - 2./D) \
         * trace(Q**2) \
         - (trace(Q)**2.) / D )
      
  b = min(1.,a)
        
  # regularized local covariance matrix for grid point 
  Qr = (1.-b)*Q + b*trace(Q)*identity(int(D)) / D
  
  return Qr

def voronoi(X,Y):
  "assign sample points to closest voronoi"
  
  n = X.shape[0]
  
  # create array associating sample points to grid points
  vlist = empty(n)
  for i,x in enumerate(X):
    dx = linalg.norm(Y-x,axis=1)
    # get index of closest distance
    vlist[i] = argmin(dx)
    if (i+1) % 10000 == 0: print " %d/%d" % (i+1,n)
  vlist = array(map(int,vlist))
  return vlist
  
def bw(X,Y,ntarget,toFile=True,doOracle=False):
  "bandwidth estimation for kernel density estimation"
  Q = cov(X.T)
  
  # check if bandwidth matrices can be stored
  if toFile:
    if not os.path.exists('bw'):
      os.makedirs('bw')
  
  N = X.shape[0]
  G = Y.shape[0]
  if len(X.shape) < 2: D = 1
  else: D = X.shape[1]
  ev = real(linalg.eigvals(Q))
  tune = max(ev)
  
  # global dimensionality
  pk = ev/sum(ev)
  pk = pk*log(pk)
  # we assume that 0*log(0) is zero
  # thus we need to check for nan values 
  # and set pk to zero for that value
  # since log(x) for x <= 0 is nan
  pk[isnan(pk)] = 0.
  d = exp(-sum(pk))
  print "  global dimensionality: %s" % d

  # set initial sigma
  sigma2 = ones(G)*tune
  
  logdetH = zeros(G)
  if not toFile: Hinv = zeros([D,D,G])
  
  for i,xi in enumerate(Y):
    # estmate localization
    w = exp(-.5*linalg.norm(X-xi,axis=1)**2/sigma2[i])
    nlocal = sum(w)
     
    # approach quickly target by using multiples of tune
    if nlocal < ntarget:
      while nlocal < ntarget:
        sigma2[i] = sigma2[i]+tune
        w = exp(-.5*linalg.norm(X-xi,axis=1)**2/sigma2[i])
        nlocal = sum(w)
     
    # fine tuning of localization approach optimal value using bisectioning
    j = 1
    while True:
      if nlocal >= ntarget: 
        sigma2[i] = sigma2[i]-tune/2.**j
      else:
        sigma2[i] = sigma2[i]+tune/2.**j
      w = exp(-.5*linalg.norm(X-xi,axis=1)**2/sigma2[i])
      nlocal = sum(w)
      if round(nlocal) == ntarget: break    
      j = j+1  
        
    # weighted covariance matrix
    Ql = cov(X.T,aweights=w)
    # estimate local dimensionality
    ev = real(linalg.eigvals(Ql))
    # global dimensionality
    pk = ev/sum(ev)
    pk = pk*log(pk)
    # we assume that 0*log(0) is zero
    # thus we need to check for nan values 
    # and set pk to zero for that value
    # since log(x) for x <= 0 is nan
    pk[isnan(pk)] = 0.
    d = exp(-sum(pk))
    
    # oracle shrinkage
    if doOracle: Ql = oracle(Ql,X.shape[0])

    # normal reference rule
    H = (4. / ( nlocal * (d+2.) ) )**( 2. / (d+4.) ) * Ql
    # logarithm of bandwidth determinant    
    sign, ldH = linalg.slogdet(H)
    logdetH[i] = ldH
    
    # store inverse bandwidth matrix
    if toFile:  save('bw/Hinv_%03d' % i,linalg.inv(H))
    else:  Hinv[:,:,i] = linalg.inv(H)

    if (i+1) % 10 == 0: print " %d/%d" % (i+1,G)

  if toFile: save('bw/logdetH',logdetH)
  else: return Hinv,logdetH
  
def logkde(X,Y,vlist,toFile=True,Hinv=None,logdetH=None):
  "estimate ln(kde) using log-sum-exp formula from numerical recipies"
   
  N = X.shape[0] # number of samples 
  if len(X.shape) < 2: d = 1
  else: d = X.shape[1]
  G = Y.shape[0]
  
  if toFile: logdetH = load('bw/logdetH.npy')
  
  kde = zeros(G)
  for i,y in enumerate(Y):
    lnKmax = 0.
    lnK = zeros(N)
    for j,x in enumerate(X):
      if i==j: continue
      # get index of voronoi associated to sample point
      v = vlist[j]
      
      # exponent of gaussian kernel at y
      dx = y-x
      if toFile: dH = matmul(dx,load('bw/Hinv_%03d.npy' % v))
      else:      dH = matmul(dx,Hinv[:,:,v])
      dHd = dot(dH,dx)
      
      # log of kernel function
      lnK[j] = -.5*(d*log(2.*pi)+logdetH[v]+dHd)
      
      # find zmax
      if lnK[j] > lnKmax: lnKmax = lnK[j]
    
    kde[i] = lnKmax + log(sum(exp(lnK-lnKmax))) - log(N)
      
    if (i+1) % 10 == 0: print " %d/%d" % (i+1,G)
    
  return kde 


X = loadtxt('traj9.cv')
ntarget = 10000
Y = X[loadtxt('smap22K.idxs',dtype=int)-1]

N = X.shape[0]
if len(X.shape) < 2: D = 1
else: D = X.shape[1]

print "voronoi association"
v = voronoi(X,Y)

print "calculating bandwidths"
Hinv,logdetH = bw(X,Y,ntarget,toFile=False,doOracle=True)
save('Hinv',Hinv)
save('logdetH',logdetH)

print "kernel density estimate"
prob = logkde(X,Y,v,toFile=False,Hinv=Hinv,logdetH=logdetH)

savetxt("lj.prob",prob.T)
  























