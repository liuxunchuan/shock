################################################################################
##                                                                            ##
##  execute this file to archive the outputs created by the MHD shock code    ##
##                                                                            ##
################################################################################

### creating the directory ###
mkdir ~/DATA_shockS/
mkdir ~/DATA_shockS/op=3.0E+00/
mkdir ~/DATA_shockS/op=3.0E+00/nH=1.0E+05cm-3/
mkdir ~/DATA_shockS/op=3.0E+00/nH=1.0E+05cm-3/Vs=2.5E+01km.s-1/

### gzip and mv output files to this directory ###
gzip output/*.out 
mv output/*.out.gz ~/DATA_shockS/op=3.0E+00/nH=1.0E+05cm-3/Vs=2.5E+01km.s-1/

### check the content of the directory ###
echo archive directory: ~/DATA_shockS/op=3.0E+00/nH=1.0E+05cm-3/Vs=2.5E+01km.s-1/
echo ---------------------------------
ls ~/DATA_shockS/op=3.0E+00/nH=1.0E+05cm-3/Vs=2.5E+01km.s-1/
