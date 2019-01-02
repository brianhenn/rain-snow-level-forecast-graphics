#!/bin/csh

cd /home/bhenn/ForecastGraphics

@ lag      = 3
echo `date --date "${lag} hours ago" '+%Y%m%d'`  >  date.txt
echo `date --date "${lag} hours ago" '+%Y'`      >> date.txt
echo `date --date "${lag} hours ago" '+%H'`      >> date.txt
echo `date --date "${lag} hours ago" '+%e'`      >> date.txt
echo `date --date "${lag} hours ago" '+%a'`      >> date.txt
echo `date --date "${lag} hours ago" '+%b'`      >> date.txt
echo `date --date "${lag} hours ago" '+%D'`      >> date.txt

set yyyy    = `date --date "${lag} hours ago" '+%Y'` 
set mm      = `date --date "${lag} hours ago" '+%m'` 
set dd      = `date --date "${lag} hours ago" '+%d'` 
set hh      = `date --date "${lag} hours ago" '+%H'`

# only use 0Z and 12Z WPC initializations with full forecast lengths
if ($hh == '06' || $hh == '18') then
echo $lag
@ lag       = $lag + 6
echo $lag
set yyyy    = `date --date "${lag} hours ago" '+%Y'`
set mm      = `date --date "${lag} hours ago" '+%m'`
set dd      = `date --date "${lag} hours ago" '+%d'`
set hh      = `date --date "${lag} hours ago" '+%H'`
endif

echo $yyyy
echo $mm
echo $dd
echo $hh

foreach fhh(006 012 018 024 030 036 042 048 054 060 066 072 078 084 090 096 102 108 114 120 126 132 138 144 150 156 162 168)
set precipfile='p06m_'$yyyy$mm$dd${hh}f$fhh'.grb'
echo $precipfile >> date.txt
ftp -n <<EOF
open ftp.wpc.ncep.noaa.gov
user anonymous bhenn@ucsd.edu
bin
prompt
cd 5km_qpf
mget $precipfile
quit
EOF
ncl_convert2nc $precipfile
rm -f p06m_*.grb
end
