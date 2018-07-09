import sys, os
import ctypes
import numpy as np
import matplotlib.pyplot as pl
from numpy.ctypeslib import ndpointer
import time


modulepath=os.path.dirname(os.path.abspath(__file__))
libfile=os.path.join(modulepath,'boris.so')

lib=ctypes.CDLL(libfile)

bwire=lib.bwire
boris_drift=lib.boris_drift
boris_kick=lib.boris_kick

arrayf64=ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")


bwire.argtypes=[
    ctypes.c_double, arrayf64, ctypes.c_int,
    arrayf64,arrayf64]

boris_drift.argtypes=[
    ctypes.c_double, arrayf64, arrayf64]

boris_kick.argtypes=[
    ctypes.c_double, arrayf64, arrayf64,
    ctypes.c_double, arrayf64, arrayf64]




clight=299792458
echarge=1.60217662e-19
class Particle(object):
    """
    E0,Ek,Pc [eV]
    q [e]
    """
    def __init__(self,u=None,gamma=None,beta=None,Ek=None,
                 vdir=[1,0,0],x=None):
        if u is not None:
          self.u=u
        elif Ek is not None:
          self.setEk(Ek,vdir)
        elif gamma is not None:
           self.setgamma(gamma,vdir)
        elif beta is not None:
           self.setbeta(beta,vdir)
        self.x=x
    def setEk(self,Ek,vdir=[1,0,0]):
        self.Ek=Ek
        gamma=Ek/self.E0+1
        beta=np.sqrt(1-1/gamma**2)
        self.setu(beta,gamma,vdir)
        return self
    def setgamma(self,gamma,vdir=[1,0,0]):
        beta=np.sqrt(1-1/gamma**2)
        self.setu(beta,gamma,vdir)
        return self
    def setbeta(self,beta,vdir=[1,0,0]):
        gamma=1/np.sqrt(1-beta**2)
        self.setu(beta,gamma,vdir)
        return self
    def setu(self,beta,gamma,vdir):
        self.beta=beta
        self.gamma=gamma
        self.u=np.array([1,vdir[0]*beta,vdir[1]*beta,vdir[2]*beta])
        self.u*=gamma
        self.E=gamma*self.E0
        self.Pc=beta*gamma*self.E0
        self.Ek=self.E-self.E0
        return self
    def __repr__(self):
        fmt="{name}(beta={self.beta},gamma={self.gamma})"
        return fmt.format(self)


class Electron(Particle):
    E0=510998.
    charge=-1

class Trajectory(object):
    def __init__(self,part,nstep,dt):
        self.x=np.zeros((nstep+1,4),dtype=float)
        self.u=np.zeros((nstep+1,4),dtype=float)
        self.x[0]=part.x
        self.u[0]=part.u
        self.q2m=part.charge/part.E0*clight**2
        self.E0=part.E0
        self.dt=dt
        self.nstep=nstep
    def track(self,bfields=[],efields=[]):
        xx=self.x[0].copy()
        uu=self.u[0].copy()
        bfield=np.zeros(3,dtype=float)
        efield=np.zeros(3,dtype=float)
        q2m=self.q2m; dt=self.dt; dt2=dt/2
        for n in range(1,self.nstep+1):
            bfield[:]=0
            efield[:]=0
            boris_drift(dt2,xx,uu)
            for src in efields:
              efield+=src.field(xx)
            for src in bfields:
              bfield+=src.field(xx)
            boris_kick(dt,xx,uu,q2m,efield,bfield)
            boris_drift(dt2,xx,uu)
            self.x[n]=xx
            self.u[n]=uu
    def plot3d(self):
        t,x,y,z=self.x.T


class BWire(object):
    def __init__(self,curr,pos):
        self.curr=curr
        self.pos=pos
    def field(self,x,bfield):
        _bfield=np.zeros(3,dtype=float)
        bwire(self.curr,self.pos,len(self.pos),x[1:4],_bfield)
        return _bfield

class BConst(object):
    def __init__(self,bx,by,bz):
        self._bfield=np.array([bx,by,bz],dtype=float)
    def field(self,x):
        return self._bfield

class EConst(object):
    def __init__(self,ex,ey,ez):
        self._efield=np.array([ex,ey,ez],dtype=float)
    def field(self,x):
        return self._efield

class FieldMap(object):
    def __init__(self,pdata,fdata):
        self.pdata=pdata
        self.fdata=fdata

def _rot(angle,i,j):
    angle=angle/180*pi
    c=cos(angle); s=sin(angle)
    mat=np.zeros((3,3),dtype=float)
    mat[i,i]=c; mat[i,j]=-s;
    mat[j,i]=s; mat[j,j]= c;
    return mat

class Move(object):
    def __init__(self,source,dx,dy,dz,xrot,yrot,zrot):
        self.source=source
        rmat=np.dot(_rot(yrot,0,2),_rot(xrot,1,2))
        rmat=np.dot(_rot(zrot,0,1),rmat)
        self.rmat=rmat
        self.svec=np.array([dx,dy,dz])
    def field(self,x):
        v=self.source(self.source.field(x))
        return np.dot(self.rmat,v+self.svec)

class Helmoltz(BWire):
   def __init__(self,cur=1,zz=0.04, rr=0.05):
       phi=np.linspace(0,2*pi,1000);
       pw1=c_[rr*cos(phi),rr*sin(phi),0*phi+zz]
       pw2=c_[rr*cos(phi),rr*sin(phi),0*phi-zz]
       self.cur=cur
       self.pos=np.vstack(pw1,pw2)



