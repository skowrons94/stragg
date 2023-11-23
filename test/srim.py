import stragg

# Get the Carbon SRIM table
table = stragg.SRIMTable( 6 )
nentries = table.getTableSize( )

# Save it to a file
f = open( "srim_table.dat", "w" )
for i in range( nentries ):
    f.write( "{}\t{}\n".format( table.getEnergy( i ), table.getStoppingPower( i ) ) )
f.close( )