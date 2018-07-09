import numpy as np
import matplotlib.pyplot as pl
import time


import pyboris

p=pyboris.Electron(beta=0.2,x=[0,0,0,0],vdir=[1,0,0])

dt=1e-9;nstep=1000

tr=pyboris.Trajectory(p,nstep,dt)

tr.track(bfields=[pyboris.BConst(0,0,1e-3)])

t,x,y,z=tr.x.T

pl.plot(x,y,'-')

pl.show()


