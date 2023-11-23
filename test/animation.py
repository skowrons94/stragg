import stragg

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation

layer = stragg.Layer( 6 )

EM = 1791
GV = 0
VV = 0
DEM = 0
Xi = 0
K = 0

DEML = layer.GetDEML(EM)
GVL  = layer.GetGVL(EM)
VVL  = layer.GetVVL(EM)
KL   = layer.GetK(EM)
Beta = layer.GetBeta(EM)
XiL  = layer.GetXi(EM)

landau = stragg.LandauFunction( )
gauss  = stragg.GaussFunction( )
vavilovAiry = stragg.VavilovAiryFunction( )
vavilovMoyal = stragg.VavilovMoyalFunction( )
vavilovEdgew = stragg.VavilovEdgeworthFunction( )

layer.setThicknessStep( 1 )

def init():
    line.set_data([], [])
    return line,

def animate(i):
    global EM, DEM, K, GV, VV, Xi
    DEML = layer.GetDEML(EM)
    GVL  = layer.GetGVL(EM)
    VVL  = layer.GetVVL(EM)
    KL   = layer.GetK(EM)
    Beta = layer.GetBeta(EM)
    XiL  = layer.GetXi(EM)

    DEM = DEM + DEML
    EM = EM - DEML
    K = K + KL
    GV = np.sqrt(GV*GV + GVL*GVL)
    VV = np.sqrt(VV*VV + VVL*VVL)
    Xi = Xi + XiL

    landau.SetLandauStep( Xi, Beta, K, DEM, 500, False )
    gauss.SetGaussStep( DEM, VV, 1000, False )
    vavilovAiry.SetAiryStep( Xi, Beta, K, DEM, 100, False )
    vavilovMoyal.SetMoyalStep( Xi, Beta, K, DEM, 100, False )
    vavilovEdgew.SetEdgeworthStep( Xi, Beta, K, DEM, 100, False )

    x = []
    land = []
    gaus = []
    vavA = []
    vavM = []
    vavE = []
    dist = []
    for k in range( 1000 ):
        x.append( 0.01*k )
        land.append( landau.GetValue( 0.01*k ) )
        gaus.append( gauss.GetValue( 0.01*k ) )
        vavA.append( vavilovAiry.GetValue( 0.01*k ) )
        vavM.append( vavilovMoyal.GetValue( 0.01*k ) )
        vavE.append( vavilovEdgew.GetValue( 0.01*k ) )
        if( K < 0.02 ):
            dist.append( landau.GetValue( 0.01*k ) )
        elif( K >= 0.02 and K < 0.30 ):
            dist.append( vavilovMoyal.GetValue( 0.01*k ) )
        elif( K >= 0.30 and k < 22.00 ):
            dist.append( vavilovAiry.GetValue( 0.01*k ) )
        elif( K>=22.00 and K < 22.00 ):
            dist.append( vavilovEdgew.GetValue( 0.01*k ) )
        else:
            dist.append( gauss.GetValue( 0.01*k ) )

    dist /= np.sum( dist )
    line.set_data(x, dist)
    ax.set_ylim(0, 1.2*max(dist))
    return line,

fig = plt.figure()
ax = plt.axes(xlim=(0, 10), ylim=(0, 1))

ax.set_xlabel( 'Energy Loss (MeV)' )
ax.set_ylabel( 'Probability' )

line, = ax.plot([], [])
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=100, interval=20, blit=True)

plt.show( )