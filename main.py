######### MAIN FILE ###########################################################################
### Integration structures for Sn photonics (pick-n-placing)
### Updated on: 10-12-2019
### Kevin C. Chen
### Notes:
###############################################################################################

import gdspy
import numpy as np
import sys

from wg_array import *


### Define parameters
unit_ys = [2.5,3,3.5]
W_inc = 0.01
W_SnV = 0.270+0.03 # for SnV
Ws_SnV = np.arange(W_SnV-2*W_inc,W_SnV+(2+1)*W_inc,W_inc)
a = 0.201783
r = 0.289494*a-0.01 # number
r_inc = 0.005
rs = np.arange(r-2*r_inc,r+(2+1)*r_inc,r_inc)
taper_L = 5
span = 4
wg_L = 15
sep_x = 1
mindist = 0.0590476
sigma = 3.95895
N_holes_wg = 5 # actual number of holes = 2*N_holes+1
N_holes_cav = 20
off = 1 # x shift inwards

N_wgs=4 # actual number of waveguides = 2*N_wgs

### Assemble
total=gdspy.Cell('main')
x0=0
y0=0

x=x0
y=y0


# ### TEST ###
# length = 2*(taper_L+sep_x+span)+wg_L
# int_struct=build_int_struct(N_wgs,unit_y,W_SnV,a,r,N_holes,taper_L,span,wg_L,sep_x,off,name='')
# total.add(gdspy.CellReference(int_struct,(x,y)))
# total.add(gdspy.Text('W='+str(int(W_SnV*1000))+' d='+str(int(r*2*1000)),2,(x-length/2-1,y-3*unit_y),layer=1,angle=np.pi/2))
# ############

# Sweep of integration waveguide structures
N=len(Ws_SnV)
for ind in range(0,N):
    x = x0
    for r in rs:
        for unit_y in unit_ys:
            W_SnV = Ws_SnV[ind]

            length = 2*(taper_L+sep_x+span)+wg_L
            label = 'Sep '+str(unit_y)+' W='+str(int(W_SnV*1000))+' d='+str(int(r*2*1000))
            wg_struct=build_wg_struct(N_wgs,unit_y,W_SnV,a,r,N_holes_wg,taper_L,span,wg_L,sep_x,off,label)
            total.add(gdspy.CellReference(wg_struct,(x,y)))
            total.add(gdspy.Text(label,2,(x-length/2.0-1,y-6.5*unit_y),layer=1,angle=np.pi/2.0))

            y = y+(2*N_wgs+5)*unit_y+30

        x = x+length+30
        y = y0

    if ind < N/2-1:
        x0 = x+10
    elif ind == N/2-1:
        x0 = 0
        y0 = y+200
        y = y0
    else: # if ind > N/2-1
        x0 = x+10

x0=0
y0=y0+200

x=x0
y=y0

# Sweep of integration cavity structures
N=len(Ws_SnV)
for ind in range(0,N):
    x = x0
    for r in rs:
        for unit_y in unit_ys:
            W_SnV = Ws_SnV[ind]

            length = 2*(taper_L+sep_x+span)+wg_L
            label = 'cSep '+str(unit_y)+' W='+str(int(W_SnV*1000))+' d='+str(int(r*2*1000))
            cav_struct=build_cav_struct(N_wgs,unit_y,W_SnV,a,r,N_holes_cav,taper_L,span,wg_L,sep_x,off,mindist,sigma,label)
            total.add(gdspy.CellReference(cav_struct,(x,y)))
            total.add(gdspy.Text(label,2,(x-length/2.0-1,y-6.5*unit_y),layer=1,angle=np.pi/2.0))

            y = y+(2*N_wgs+5)*unit_y+30

        x = x+length+30
        y = y0

    if ind < N/2-1:
        x0 = x+10
    elif ind == N/2-1:
        x0 = 0
        y0 = y+200
        y = y0
    else: # if ind > N/2-1
        x0 = x+10

gdspy.write_gds('SnV_photonics_test_v2.gds')
