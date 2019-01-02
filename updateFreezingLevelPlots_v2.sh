#!/bin/bash

# Script to download, process, and update freezing level forecast plots for the GEFS visualizations shown on the CW3E/DWR/AR decision support website
# Brian Henn, CW3E/SIO/UCSD, December 2017
# bhenn@ucsd.edu
# v2 created in January 2017 to use WPC 5km precipitation instead of GFS

source .bashrc
cd ~/ForecastGraphics/

# download grib files of freezing levels from GEFS/NOMADS
./fzl.sc > fzlLog.txt

# download grib files of precipitation from GEFS/NOMADS
#./qpf.sc > qpfLog.txt

# download grib files from WPC FTP and convert to netCDFs
rm *.nc
./wpc_qpf.sc > qpfLog.txt

# run script to covert grib files to netCDFs
ncl ./convertgrb2nc.ncl > convert1log.txt
#ncl ./convertPgrb2nc.ncl > convert2log.txt

# now call matlab script to create figures
rm ./outputFiles/*.png 
cd ./matlab/
matlab -nodisplay -nosplash -r plotGEFSfreezingLevelsv4_skyriver > matlaberr.log $2> matlabout.log
#mv ./outputFiles/*current


# sync output files to CW3E web server
cd ..
rsync -r --delete ./outputFiles/ cw3e@cw3e-web.ucsd.edu:/var/www/cw3e/htdocs/images/gefs/freezingLevelImages/ 
