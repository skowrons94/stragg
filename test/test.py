import stragg

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

DEML = layer.GetDEML(EM)
GVL  = layer.GetGVL(EM)
VVL  = layer.GetVVL(EM)
KL   = layer.GetK(EM)
Beta = layer.GetBeta(EM)
XiL  = layer.GetXi(EM)

landau.SetLandauStep( Xi, Beta, K, DEM, 500, False )
gauss.SetGaussStep( DEM, VV, 1000, False )
vavilovAiry.SetAiryStep( Xi, Beta, K, DEM, 100, False )
vavilovMoyal.SetMoyalStep( Xi, Beta, K, DEM, 100, False )
vavilovEdgew.SetEdgeworthStep( Xi, Beta, K, DEM, 100, False )