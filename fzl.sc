#!/bin/csh

#cd ~/Models/ARPortal/Requests/fzl/

#set yyyy = 2017
#set mm   = 02
#set dd   = 14
#set hh   = 12

#set lag = 10
set lag = 0
echo `date --date "${lag} hours ago" '+%Y%m%d'`  >  date.txt
echo `date --date "${lag} hours ago" '+%Y'`      >> date.txt
echo `date --date "${lag} hours ago" '+%H'`      >> date.txt
echo `date --date "${lag} hours ago" '+%e'`      >> date.txt
echo `date --date "${lag} hours ago" '+%a'`      >> date.txt
echo `date --date "${lag} hours ago" '+%b'`      >> date.txt
echo `date --date "${lag} hours ago" '+%D'`      >> date.txt

set date    = `date --date "${lag} hours ago" '+%Y%m%d'`
set hh      = `date --date "${lag} hours ago" '+%H'`

#echo $date
#echo $hh

set get_data = 1

if ($get_data == 1) then 

'rm' out.grb
foreach ens(c00 p01 p02 p03 p04 p05 p06 p07 p08 p09 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20)
foreach fhh(000 003 006 009 012 015 018 021 024 027 030 033 036 039 042 045 048 051 054 057 060 063 066 069 072 075 078 081 084 087 090 093 096 099 102 105 108 111 114 117 120 123 126 129 132 135 138 141 144 147 150 153 156 159 162 165 168 171 174 177 180 183 186 189 192)
#foreach fhh(000 003 006 009 012)
set idx = "http://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/gefs."${date}"/"${hh}"/pgrb2bp5/ge"${ens}".t"${hh}"z.pgrb2b.0p50.f"${fhh}".idx"
echo $idx
set url = "http://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/gefs."${date}"/"${hh}"/pgrb2bp5/ge"${ens}".t"${hh}"z.pgrb2b.0p50.f"${fhh}
./get_inv.pl ${idx} | grep ":HGT:0C" | ./get_grib.pl ${url} out_${fhh}_${ens}.grb
end 
end
cat *.grb > out.grb
'mv' out.grb out.grb2
'rm' *.grb 
'mv' out.grb2 out.grb

endif

#foreach thresh(7000 7500 8000 8500 9000 9500 10000)

#echo $thresh > thresh.txt

#ncl fzl_prob.ncl

#end

#set path = (/usr/local/bin $path)
#setenv NCARG_ROOT /usr/local

#ncl landfall_plumes_7day.ncl

exit(0)
