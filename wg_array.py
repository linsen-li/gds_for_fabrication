import gdspy
import numpy as np
import sys


##############################
# Builds waveguide unit cell
##############################
def wg(unit_y,W,a,r,N,taper_L,span,wg_L,sep_x,off,name=''):
    ##########
    # unit_y:       length per unit array cell in y
    # W:            waveguide width
    # a:            periodicty of Bragg mirrors at the waveguide's center
    # r:            radius of circles
    # N:            number of Bragg mirrors
    # taper_L:      length over which linear tapering occurs
    # span:         x span of Gaussian bulge at the interconnection
    # wg_L:         length of the straight waveguide section
    # sep_x:        separation from the bulk
    ##########

    cell_wg = gdspy.Cell('cell_wg_'+name)

    ## Add big trench for subtraction later
    length=2*(sep_x+taper_L+span)+wg_L
    trench = gdspy.Rectangle((-length/2.0,unit_y/2.0),(length/2.0,-unit_y/2.0),layer=1)

    ## Waveguide
    straight_wg = gdspy.Rectangle((-wg_L/2.0-span,-W/2),(wg_L/2.0+span,W/2),layer=1)

    ## Connecting vertical beams
    width=1.2*W # width of the support beams
    pos_x = -span/2.0-wg_L/2.0+off
    points = [(pos_x-width/2,-unit_y/2.0),(pos_x-width/2.0,-width-W),(pos_x,0),\
        (pos_x-width/2,width+W),(pos_x-width/2,unit_y/2.0),(pos_x+width/2,unit_y/2.0),\
        (pos_x+width/2,width+W),(pos_x,0),(pos_x+width/2.0,-width-W),\
        (pos_x+width/2,-unit_y/2.0)]
    beam_L = gdspy.Polygon(points,layer=1)
    pos_x = span/2.0+wg_L/2.0-off
    points = [(pos_x-width/2,-unit_y/2.0),(pos_x-width/2.0,-width-W),(pos_x,0),\
        (pos_x-width/2,width+W),(pos_x-width/2,unit_y/2.0),(pos_x+width/2,unit_y/2.0),\
        (pos_x+width/2,width+W),(pos_x,0),(pos_x+width/2.0,-width-W),\
        (pos_x+width/2,-unit_y/2.0)]
    beam_R = gdspy.Polygon(points,layer=1)

    ## Gaussian-shaped intersecitons
    amplitude = 0.1
    sd = 0.6
    p = np.linspace(-span/2.0,span/2.0,199)
    q = amplitude*np.exp(-p**2/2/sd/sd)

    x_offset = wg_L/2.0+span/2.0-off
    gaussLT = gdspy.Polygon(list(zip(p-x_offset,q+W/2.0)))
    gaussLB = gdspy.Polygon(list(zip(p-x_offset,-q-W/2.0)))
    gaussRT = gdspy.Polygon(list(zip(p+x_offset,q+W/2.0)))
    gaussRB = gdspy.Polygon(list(zip(p+x_offset,-q-W/2.0)))

    ## Add tapers
    points = [(-wg_L/2.0-span,W/2.0),(-wg_L/2.0-span-taper_L,0),(-wg_L/2.0-span,-W/2.0)]
    taperL = gdspy.Polygon(points,layer=1)
    points = [(wg_L/2.0+span,W/2.0),(wg_L/2.0+span+taper_L,0),(wg_L/2.0+span,-W/2.0)]
    taperR = gdspy.Polygon(points,layer=1)

    ## Add holes
    pos = 0
    hole = gdspy.Round((pos,0),r,layer=1)
    cell_wg.add(hole)
    for n in range(0,N):
        pos = pos+a
        hole = gdspy.Round((pos,0),r,layer=1)
        cell_wg.add(hole)
        hole = gdspy.Round((-pos,0),r,layer=1)
        cell_wg.add(hole)

    ## Boolean operations
    beams = gdspy.fast_boolean(beam_L,beam_R,'or')
    gaussL = gdspy.fast_boolean(gaussLT,gaussLB,'or')
    gaussR = gdspy.fast_boolean(gaussRT,gaussRB,'or')
    gauss = gdspy.fast_boolean(gaussL,gaussR,'or')
    tapers = gdspy.fast_boolean(taperL,taperR,'or')

    sum1a = gdspy.fast_boolean(beams,gauss,'or')
    sum1b = gdspy.fast_boolean(straight_wg,tapers,'or')
    sum1 = gdspy.fast_boolean(sum1a,sum1b,'or')

    cell_wg.add(gdspy.fast_boolean(trench,sum1,'not'))

    return cell_wg


##############################
# Builds waveguide array
##############################
def build_wg_array(N_wgs,unit_y,W,a,r,N,taper_L,span,wg_L,sep_x,off,name=''):

    wg_array = gdspy.Cell('wg_array_'+name)

    x0=0
    y0=unit_y/2.0

    # Build array
    wg_cell = wg(unit_y,W,a,r,N,taper_L,span,wg_L,sep_x,off,name)

    for n in range(0,N_wgs):
        wg_array.add(gdspy.CellReference(wg_cell,(x0,y0)))
        wg_array.add(gdspy.CellReference(wg_cell,(x0,-y0)))
        y0 = y0+unit_y

    return wg_array


#################################################
# Build support frame structures for waveguides
#################################################
def makeFrame_wg(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,name=''):
    # Notes: bottom frame (origin (0,0) set at the center of the bottom beam in the frame)

    frame = gdspy.Cell('frame_wg_'+name)

    ############### Bottom

    ### define trench for later Boolean subtraction
    length = 2*(sep_x+taper_L+span)+wg_L
    trench_large = gdspy.Rectangle((-length/2.0,0),(length/2.0,-unit_y/2.0),layer=1)

    length = 2*span+wg_L
    trench_medium = gdspy.Rectangle((-length/2.0,-unit_y/2.0),(length/2.0,-1.5*unit_y),layer=1)

    posx_holes = 0
    for n in range(0,N_holes):
        posx_holes = posx_holes+a
    length = posx_holes+unit_y*2
    trench_small = gdspy.Rectangle((-length/2.0,-1.5*unit_y),(length/2.0,-2.5*unit_y),layer=1)

    trench1 = gdspy.fast_boolean(trench_large,trench_medium,'or')
    trench = gdspy.fast_boolean(trench_small,trench1,'or')

    #### define bones
    # L braces
    width = 2*W # width of the L brace
    posx = -(wg_L/2.0+span/2.0-off-W/2.0)
    points = [(posx,0),(posx-width,0),(posx-width,-unit_y/2.0-width),\
        (posx+unit_y/2.0+width,-unit_y/2.0-width),(posx+unit_y/2.0+width,-unit_y/2.0),\
            (posx,-unit_y/2.0)]
    L_left = gdspy.Polygon(points,layer=1)

    beamwidth = 0.75*W # width of the support beam for the brace
    points = [(posx,-beamwidth),(posx,-2*beamwidth),\
        (posx+unit_y/2.0+width-2*beamwidth,-unit_y/2.0),(posx+unit_y/2.0+width-beamwidth,-unit_y/2.0)]
    L_left_supp = gdspy.Polygon(points,layer=1)

    posx = wg_L/2.0+span/2.0-off-W/2.0
    points = [(posx,0),(posx+width,0),(posx+width,-unit_y/2.0-width),\
        (posx-unit_y/2.0-width,-unit_y/2.0-width),(posx-unit_y/2.0-width,-unit_y/2.0),\
            (posx,-unit_y/2.0)]
    L_right = gdspy.Polygon(points,layer=1)

    points = [(posx,-beamwidth),(posx,-2*beamwidth),\
        (posx-unit_y/2.0-width+2*beamwidth,-unit_y/2.0),(posx-unit_y/2.0-width+beamwidth,-unit_y/2.0)]
    L_right_supp = gdspy.Polygon(points,layer=1)

    brace_left = gdspy.fast_boolean(L_left,L_left_supp,'or')
    brace_right = gdspy.fast_boolean(L_right,L_right_supp,'or')
    braces = gdspy.fast_boolean(brace_left,brace_right,'or')

    # horizontal and vertical beams
    beam_h_1a = gdspy.Rectangle((posx-unit_y/2.0-width,-unit_y/2.0),(-posx+unit_y/2.0+width,-unit_y/2.0-W),layer=1)

    width = 0.9*W
    posx = wg_L/2.0-2*off
    beam_v_1a = gdspy.Rectangle((posx,-unit_y/2.0-W),(posx-width,-1.5*unit_y),layer=1)
    beam_v_1b = gdspy.Rectangle((-posx,-unit_y/2.0-W),(-posx+width,-1.5*unit_y),layer=1)
    beam_v_1 = gdspy.fast_boolean(beam_v_1a,beam_v_1b,'or')

    beam_h_1b = gdspy.Rectangle((posx_holes,-1.5*unit_y),(-posx_holes,-1.5*unit_y+W),layer=1)
    beam_h_1 = gdspy.fast_boolean(beam_h_1a,beam_h_1b,'or')

    beam_v_2a = gdspy.Rectangle((posx_holes,-unit_y/2.0-W),(posx_holes-W,-1.5*unit_y),layer=1)
    beam_v_2b = gdspy.Rectangle((-posx_holes,-unit_y/2.0-W),(-posx_holes+W,-1.5*unit_y),layer=1)
    beam_v_2 = gdspy.fast_boolean(beam_v_2a,beam_v_2b,'or')

    beam_v_3a = gdspy.Rectangle((posx_holes,-1.5*unit_y),(posx_holes-width,-2.5*unit_y),layer=1)
    beam_v_3b = gdspy.Rectangle((-posx_holes,-1.5*unit_y),(-posx_holes+width,-2.5*unit_y),layer=1)
    beam_v_3 = gdspy.fast_boolean(beam_v_3a,beam_v_3b,'or')

    beam_v_temp = gdspy.fast_boolean(beam_v_1,beam_v_2,'or')
    beam_v = gdspy.fast_boolean(beam_v_temp,beam_v_3,'or')
    beams = gdspy.fast_boolean(beam_v,beam_h_1,'or')

    bones_bot = gdspy.fast_boolean(braces,beams,'or')

    frame.add(gdspy.fast_boolean(trench,bones_bot,'not'))

    ############### Top

    ### define trenches
    y0=8*unit_y
    length = 2*(sep_x+taper_L+span)+wg_L
    trench_large = gdspy.Rectangle((-length/2.0,y0),(length/2.0,y0+unit_y/2.0),layer=1)

    length = 2*span+wg_L
    trench_medium = gdspy.Rectangle((-length/2.0,y0+unit_y/2.0),(length/2.0,y0+1.5*unit_y),layer=1)

    length = posx_holes+unit_y*2
    trench_small = gdspy.Rectangle((-length/2.0,y0+1.5*unit_y),(length/2.0,y0+2.5*unit_y),layer=1)

    trench1 = gdspy.fast_boolean(trench_large,trench_medium,'or')
    trench = gdspy.fast_boolean(trench_small,trench1,'or')

    #### define bones
    # L braces
    width = 2*W # width of the L brace
    posx = -(wg_L/2.0+span/2.0-off-W/2.0)
    points = [(posx,y0),(posx-width,y0),(posx-width,y0+unit_y/2.0+width),\
        (posx+unit_y/2.0+width,y0+unit_y/2.0+width),(posx+unit_y/2.0+width,y0+unit_y/2.0),\
            (posx,y0+unit_y/2.0)]
    L_left = gdspy.Polygon(points,layer=1)

    beamwidth = 0.75*W # width of the support beam for the brace
    points = [(posx,y0+beamwidth),(posx,y0+2*beamwidth),\
        (posx+unit_y/2.0+width-2*beamwidth,y0+unit_y/2.0),(posx+unit_y/2.0+width-beamwidth,y0+unit_y/2.0)]
    L_left_supp = gdspy.Polygon(points,layer=1)

    posx = wg_L/2.0+span/2.0-off-W/2.0
    points = [(posx,y0),(posx+width,y0),(posx+width,y0+unit_y/2.0+width),\
        (posx-unit_y/2.0-width,y0+unit_y/2.0+width),(posx-unit_y/2.0-width,y0+unit_y/2.0),\
            (posx,y0+unit_y/2.0)]
    L_right = gdspy.Polygon(points,layer=1)

    points = [(posx,y0+beamwidth),(posx,y0+2*beamwidth),\
        (posx-unit_y/2.0-width+2*beamwidth,y0+unit_y/2.0),(posx-unit_y/2.0-width+beamwidth,y0+unit_y/2.0)]
    L_right_supp = gdspy.Polygon(points,layer=1)

    brace_left = gdspy.fast_boolean(L_left,L_left_supp,'or')
    brace_right = gdspy.fast_boolean(L_right,L_right_supp,'or')
    braces = gdspy.fast_boolean(brace_left,brace_right,'or')

    # horizontal and vertical beams
    beam_h_1a = gdspy.Rectangle((posx-unit_y/2.0-width,y0+unit_y/2.0),(-posx+unit_y/2.0+width,y0+unit_y/2.0+W),layer=1)

    width = 0.9*W
    posx = wg_L/2.0-2*off
    beam_v_1a = gdspy.Rectangle((posx,y0+unit_y/2.0+W),(posx-width,y0+1.5*unit_y),layer=1)
    beam_v_1b = gdspy.Rectangle((-posx,y0+unit_y/2.0+W),(-posx+width,y0+1.5*unit_y),layer=1)
    beam_v_1 = gdspy.fast_boolean(beam_v_1a,beam_v_1b,'or')

    beam_h_1b = gdspy.Rectangle((posx_holes,y0+1.5*unit_y),(-posx_holes,y0+1.5*unit_y+W),layer=1)
    beam_h_1 = gdspy.fast_boolean(beam_h_1a,beam_h_1b,'or')

    points = [(-posx_holes+W,y0+unit_y/2.0+W),(-posx_holes+W+width,y0+unit_y/2.0+W),\
        (posx_holes-W,y0+1.5*unit_y-width),(posx_holes-W,y0+1.5*unit_y)]
    beam_supp = gdspy.Polygon(points,layer=1)
    beam_h = gdspy.fast_boolean(beam_h_1,beam_supp,'or')

    beam_v_2a = gdspy.Rectangle((posx_holes,y0+unit_y/2.0+W),(posx_holes-W,y0+1.5*unit_y),layer=1)
    beam_v_2b = gdspy.Rectangle((-posx_holes,y0+unit_y/2.0+W),(-posx_holes+W,y0+1.5*unit_y),layer=1)
    beam_v_2 = gdspy.fast_boolean(beam_v_2a,beam_v_2b,'or')

    beam_v_3a = gdspy.Rectangle((posx_holes,y0+1.5*unit_y),(posx_holes-width,y0+2.5*unit_y),layer=1)
    beam_v_3b = gdspy.Rectangle((-posx_holes,y0+1.5*unit_y),(-posx_holes+width,y0+2.5*unit_y),layer=1)
    beam_v_3 = gdspy.fast_boolean(beam_v_3a,beam_v_3b,'or')

    beam_v_temp = gdspy.fast_boolean(beam_v_1,beam_v_2,'or')
    beam_v = gdspy.fast_boolean(beam_v_temp,beam_v_3,'or')
    beams = gdspy.fast_boolean(beam_v,beam_h,'or')

    bones_top = gdspy.fast_boolean(braces,beams,'or')

    frame.add(gdspy.fast_boolean(trench,bones_top,'not'))

    return frame


##############################
# Combine frame and waveguide array
##############################
def build_wg_struct(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,name=''):

    wg_struct = gdspy.Cell('integration_wg_'+name)

    wg_array = build_wg_array(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,name)
    frame = makeFrame_wg(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,name)

    x=0
    y=0
    wg_struct.add(gdspy.CellReference(wg_array,(x,y)))

    y=-N_wgs*unit_y
    wg_struct.add(gdspy.CellReference(frame,(x,y)))

    return wg_struct


##############################
# Builds cavity unit cell
##############################
def cav(unit_y,W,a,r,N,taper_L,span,wg_L,sep_x,off,mindist,sigma,name=''):
    ##########
    # unit_y:       length per unit array cell in y
    # W:            waveguide width
    # a:            periodicty of Bragg mirrors at the waveguide's center
    # r:            radius of circles
    # N:            number of Bragg mirrors
    # taper_L:      length over which linear tapering occurs
    # span:         x span of Gaussian bulge at the interconnection
    # wg_L:         length of the straight waveguide section
    # sep_x:        separation from the bulk
    ##########

    cell_cav = gdspy.Cell('cell_cav_'+name)

    ## Add big trench for subtraction later
    length=2*(sep_x+taper_L+span)+wg_L
    trench = gdspy.Rectangle((-length/2.0,unit_y/2.0),(length/2.0,-unit_y/2.0),layer=1)

    ## Waveguide
    straight_wg = gdspy.Rectangle((-wg_L/2.0-span,-W/2),(wg_L/2.0+span,W/2),layer=1)

    ## Connecting vertical beams
    width=1.0*W # width of the support beams
    pos_x = -span/2.0-wg_L/2.0+off
    points = [(pos_x-width/2,-unit_y/2.0),(pos_x-width/2.0,-width-W),(pos_x,0),\
        (pos_x-width/2,width+W),(pos_x-width/2,unit_y/2.0),(pos_x+width/2,unit_y/2.0),\
        (pos_x+width/2,width+W),(pos_x,0),(pos_x+width/2.0,-width-W),\
        (pos_x+width/2,-unit_y/2.0)]
    beam_L = gdspy.Polygon(points,layer=1)
    pos_x = span/2.0+wg_L/2.0-off
    points = [(pos_x-width/2,-unit_y/2.0),(pos_x-width/2.0,-width-W),(pos_x,0),\
        (pos_x-width/2,width+W),(pos_x-width/2,unit_y/2.0),(pos_x+width/2,unit_y/2.0),\
        (pos_x+width/2,width+W),(pos_x,0),(pos_x+width/2.0,-width-W),\
        (pos_x+width/2,-unit_y/2.0)]
    beam_R = gdspy.Polygon(points,layer=1)

    ## Gaussian-shaped intersecitons
    amplitude = 0.1
    sd = 0.6
    p = np.linspace(-span/2.0,span/2.0,199)
    q = amplitude*np.exp(-p**2/2/sd/sd)

    x_offset = wg_L/2.0+span/2.0-off
    gaussLT = gdspy.Polygon(list(zip(p-x_offset,q+W/2.0)))
    gaussLB = gdspy.Polygon(list(zip(p-x_offset,-q-W/2.0)))
    gaussRT = gdspy.Polygon(list(zip(p+x_offset,q+W/2.0)))
    gaussRB = gdspy.Polygon(list(zip(p+x_offset,-q-W/2.0)))

    ## Add tapers
    points = [(-wg_L/2.0-span,W/2.0),(-wg_L/2.0-span-taper_L,0),(-wg_L/2.0-span,-W/2.0)]
    taperL = gdspy.Polygon(points,layer=1)
    points = [(wg_L/2.0+span,W/2.0),(wg_L/2.0+span+taper_L,0),(wg_L/2.0+span,-W/2.0)]
    taperR = gdspy.Polygon(points,layer=1)

    ## Add holes
    a0=mindist+2*r
    pos=-a0/2

    for n in range(0,N):
        aa=a-(a-a0)*np.exp(-n**2/(2*sigma**2))
        pos=pos+aa
        hole = gdspy.Round((pos,0),r,layer=1)
        cell_cav.add(hole)
        hole = gdspy.Round((-pos,0),r,layer=1)
        cell_cav.add(hole)


    ## Boolean operations
    beams = gdspy.fast_boolean(beam_L,beam_R,'or')
    gaussL = gdspy.fast_boolean(gaussLT,gaussLB,'or')
    gaussR = gdspy.fast_boolean(gaussRT,gaussRB,'or')
    gauss = gdspy.fast_boolean(gaussL,gaussR,'or')
    tapers = gdspy.fast_boolean(taperL,taperR,'or')

    sum1a = gdspy.fast_boolean(beams,gauss,'or')
    sum1b = gdspy.fast_boolean(straight_wg,tapers,'or')
    sum1 = gdspy.fast_boolean(sum1a,sum1b,'or')

    cell_cav.add(gdspy.fast_boolean(trench,sum1,'not'))

    return cell_cav


##############################
# Builds cavity array
##############################
def build_cav_array(N_wgs,unit_y,W,a,r,N,taper_L,span,wg_L,sep_x,off,mindist,sigma,name=''):

    cav_array = gdspy.Cell('cav_array_'+name)

    x0=0
    y0=unit_y/2.0

    # Build array
    cav_cell = cav(unit_y,W,a,r,N,taper_L,span,wg_L,sep_x,off,mindist,sigma,name)

    for n in range(0,N_wgs):
        cav_array.add(gdspy.CellReference(cav_cell,(x0,y0)))
        cav_array.add(gdspy.CellReference(cav_cell,(x0,-y0)))
        y0 = y0+unit_y

    return cav_array


#################################################
# Build support frame structures for waveguides
#################################################
def makeFrame_cav(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,mindist,sigma,name=''):
    # Notes: bottom frame (origin (0,0) set at the center of the bottom beam in the frame)

    frame = gdspy.Cell('frame_cav_'+name)

    ############### Bottom

    ### define trench for later Boolean subtraction
    length = 2*(sep_x+taper_L+span)+wg_L
    trench_large = gdspy.Rectangle((-length/2.0,0),(length/2.0,-unit_y/2.0),layer=1)

    length = 2*span+wg_L
    trench_medium = gdspy.Rectangle((-length/2.0,-unit_y/2.0),(length/2.0,-1.5*unit_y),layer=1)

    posx_holes = 0
    for n in range(0,N_wgs):
        posx_holes = posx_holes+a
    length = posx_holes+unit_y*2
    trench_small = gdspy.Rectangle((-length/2.0,-1.5*unit_y),(length/2.0,-2.5*unit_y),layer=1)

    trench1 = gdspy.fast_boolean(trench_large,trench_medium,'or')
    trench = gdspy.fast_boolean(trench_small,trench1,'or')

    #### define bones
    # L braces
    width = 2*W # width of the L brace
    posx = -(wg_L/2.0+span/2.0-off-W/2.0)
    points = [(posx,0),(posx-width,0),(posx-width,-unit_y/2.0-width),\
        (posx+unit_y/2.0+width,-unit_y/2.0-width),(posx+unit_y/2.0+width,-unit_y/2.0),\
            (posx,-unit_y/2.0)]
    L_left = gdspy.Polygon(points,layer=1)

    beamwidth = 0.75*W # width of the support beam for the brace
    points = [(posx,-beamwidth),(posx,-2*beamwidth),\
        (posx+unit_y/2.0+width-2*beamwidth,-unit_y/2.0),(posx+unit_y/2.0+width-beamwidth,-unit_y/2.0)]
    L_left_supp = gdspy.Polygon(points,layer=1)

    posx = wg_L/2.0+span/2.0-off-W/2.0
    points = [(posx,0),(posx+width,0),(posx+width,-unit_y/2.0-width),\
        (posx-unit_y/2.0-width,-unit_y/2.0-width),(posx-unit_y/2.0-width,-unit_y/2.0),\
            (posx,-unit_y/2.0)]
    L_right = gdspy.Polygon(points,layer=1)

    points = [(posx,-beamwidth),(posx,-2*beamwidth),\
        (posx-unit_y/2.0-width+2*beamwidth,-unit_y/2.0),(posx-unit_y/2.0-width+beamwidth,-unit_y/2.0)]
    L_right_supp = gdspy.Polygon(points,layer=1)

    brace_left = gdspy.fast_boolean(L_left,L_left_supp,'or')
    brace_right = gdspy.fast_boolean(L_right,L_right_supp,'or')
    braces = gdspy.fast_boolean(brace_left,brace_right,'or')

    # horizontal and vertical beams
    beam_h_1a = gdspy.Rectangle((posx-unit_y/2.0-width,-unit_y/2.0),(-posx+unit_y/2.0+width,-unit_y/2.0-W),layer=1)

    width = 0.9*W
    posx = wg_L/2.0-2*off
    beam_v_1a = gdspy.Rectangle((posx,-unit_y/2.0-W),(posx-width,-1.5*unit_y),layer=1)
    beam_v_1b = gdspy.Rectangle((-posx,-unit_y/2.0-W),(-posx+width,-1.5*unit_y),layer=1)
    beam_v_1 = gdspy.fast_boolean(beam_v_1a,beam_v_1b,'or')

    beam_h_1b = gdspy.Rectangle((posx_holes,-1.5*unit_y),(-posx_holes,-1.5*unit_y+W),layer=1)
    beam_h_1 = gdspy.fast_boolean(beam_h_1a,beam_h_1b,'or')

    beam_v_2a = gdspy.Rectangle((posx_holes,-unit_y/2.0-W),(posx_holes-W,-1.5*unit_y),layer=1)
    beam_v_2b = gdspy.Rectangle((-posx_holes,-unit_y/2.0-W),(-posx_holes+W,-1.5*unit_y),layer=1)
    beam_v_2 = gdspy.fast_boolean(beam_v_2a,beam_v_2b,'or')

    beam_v_3a = gdspy.Rectangle((posx_holes,-1.5*unit_y),(posx_holes-width,-2.5*unit_y),layer=1)
    beam_v_3b = gdspy.Rectangle((-posx_holes,-1.5*unit_y),(-posx_holes+width,-2.5*unit_y),layer=1)
    beam_v_3 = gdspy.fast_boolean(beam_v_3a,beam_v_3b,'or')

    beam_v_temp = gdspy.fast_boolean(beam_v_1,beam_v_2,'or')
    beam_v = gdspy.fast_boolean(beam_v_temp,beam_v_3,'or')
    beams = gdspy.fast_boolean(beam_v,beam_h_1,'or')

    bones_bot = gdspy.fast_boolean(braces,beams,'or')

    frame.add(gdspy.fast_boolean(trench,bones_bot,'not'))

    ############### Top

    ### define trenches
    y0=8*unit_y
    length = 2*(sep_x+taper_L+span)+wg_L
    trench_large = gdspy.Rectangle((-length/2.0,y0),(length/2.0,y0+unit_y/2.0),layer=1)

    length = 2*span+wg_L
    trench_medium = gdspy.Rectangle((-length/2.0,y0+unit_y/2.0),(length/2.0,y0+1.5*unit_y),layer=1)

    length = posx_holes+unit_y*2
    trench_small = gdspy.Rectangle((-length/2.0,y0+1.5*unit_y),(length/2.0,y0+2.5*unit_y),layer=1)

    trench1 = gdspy.fast_boolean(trench_large,trench_medium,'or')
    trench = gdspy.fast_boolean(trench_small,trench1,'or')

    #### define bones
    # L braces
    width = 2*W # width of the L brace
    posx = -(wg_L/2.0+span/2.0-off-W/2.0)
    points = [(posx,y0),(posx-width,y0),(posx-width,y0+unit_y/2.0+width),\
        (posx+unit_y/2.0+width,y0+unit_y/2.0+width),(posx+unit_y/2.0+width,y0+unit_y/2.0),\
            (posx,y0+unit_y/2.0)]
    L_left = gdspy.Polygon(points,layer=1)

    beamwidth = 0.75*W # width of the support beam for the brace
    points = [(posx,y0+beamwidth),(posx,y0+2*beamwidth),\
        (posx+unit_y/2.0+width-2*beamwidth,y0+unit_y/2.0),(posx+unit_y/2.0+width-beamwidth,y0+unit_y/2.0)]
    L_left_supp = gdspy.Polygon(points,layer=1)

    posx = wg_L/2.0+span/2.0-off-W/2.0
    points = [(posx,y0),(posx+width,y0),(posx+width,y0+unit_y/2.0+width),\
        (posx-unit_y/2.0-width,y0+unit_y/2.0+width),(posx-unit_y/2.0-width,y0+unit_y/2.0),\
            (posx,y0+unit_y/2.0)]
    L_right = gdspy.Polygon(points,layer=1)

    points = [(posx,y0+beamwidth),(posx,y0+2*beamwidth),\
        (posx-unit_y/2.0-width+2*beamwidth,y0+unit_y/2.0),(posx-unit_y/2.0-width+beamwidth,y0+unit_y/2.0)]
    L_right_supp = gdspy.Polygon(points,layer=1)

    brace_left = gdspy.fast_boolean(L_left,L_left_supp,'or')
    brace_right = gdspy.fast_boolean(L_right,L_right_supp,'or')
    braces = gdspy.fast_boolean(brace_left,brace_right,'or')

    # horizontal and vertical beams
    beam_h_1a = gdspy.Rectangle((posx-unit_y/2.0-width,y0+unit_y/2.0),(-posx+unit_y/2.0+width,y0+unit_y/2.0+W),layer=1)

    width = 0.9*W
    posx = wg_L/2.0-2*off
    beam_v_1a = gdspy.Rectangle((posx,y0+unit_y/2.0+W),(posx-width,y0+1.5*unit_y),layer=1)
    beam_v_1b = gdspy.Rectangle((-posx,y0+unit_y/2.0+W),(-posx+width,y0+1.5*unit_y),layer=1)
    beam_v_1 = gdspy.fast_boolean(beam_v_1a,beam_v_1b,'or')

    beam_h_1b = gdspy.Rectangle((posx_holes,y0+1.5*unit_y),(-posx_holes,y0+1.5*unit_y+W),layer=1)
    beam_h_1 = gdspy.fast_boolean(beam_h_1a,beam_h_1b,'or')

    points = [(-posx_holes+W,y0+unit_y/2.0+W),(-posx_holes+W+width,y0+unit_y/2.0+W),\
        (posx_holes-W,y0+1.5*unit_y-width),(posx_holes-W,y0+1.5*unit_y)]
    beam_supp = gdspy.Polygon(points,layer=1)
    beam_h = gdspy.fast_boolean(beam_h_1,beam_supp,'or')

    beam_v_2a = gdspy.Rectangle((posx_holes,y0+unit_y/2.0+W),(posx_holes-W,y0+1.5*unit_y),layer=1)
    beam_v_2b = gdspy.Rectangle((-posx_holes,y0+unit_y/2.0+W),(-posx_holes+W,y0+1.5*unit_y),layer=1)
    beam_v_2 = gdspy.fast_boolean(beam_v_2a,beam_v_2b,'or')

    beam_v_3a = gdspy.Rectangle((posx_holes,y0+1.5*unit_y),(posx_holes-width,y0+2.5*unit_y),layer=1)
    beam_v_3b = gdspy.Rectangle((-posx_holes,y0+1.5*unit_y),(-posx_holes+width,y0+2.5*unit_y),layer=1)
    beam_v_3 = gdspy.fast_boolean(beam_v_3a,beam_v_3b,'or')

    beam_v_temp = gdspy.fast_boolean(beam_v_1,beam_v_2,'or')
    beam_v = gdspy.fast_boolean(beam_v_temp,beam_v_3,'or')
    beams = gdspy.fast_boolean(beam_v,beam_h,'or')

    bones_top = gdspy.fast_boolean(braces,beams,'or')

    frame.add(gdspy.fast_boolean(trench,bones_top,'not'))

    return frame


##############################
# Combine frame and cavity array
##############################
def build_cav_struct(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,mindist,sigma,name=''):

    cav_struct = gdspy.Cell('integration_cav_'+name)

    cav_array = build_cav_array(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,mindist,sigma,name)
    frame = makeFrame_cav(N_wgs,unit_y,W,a,r,N_holes,taper_L,span,wg_L,sep_x,off,mindist,sigma,name)

    x=0
    y=0
    cav_struct.add(gdspy.CellReference(cav_array,(x,y)))

    y=-N_wgs*unit_y
    cav_struct.add(gdspy.CellReference(frame,(x,y)))

    return cav_struct