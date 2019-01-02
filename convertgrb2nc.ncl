;*****************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;***************************************
begin

infile  = "/home/bhenn/ForecastGraphics/out.grb"
outfile = "./out.nc"

a = addfile(infile,"r")

lat1 = 30
lat2 = 52
lon1 = 230
lon2 = 250

z             = a->HGT_P1_L4_GLL0(:,:,{lat1:lat2},{lon1:lon2})
lat           = a->lat_0({lat1:lat2})
lon           = a->lon_0({lon1:lon2})
forecast_time = a->forecast_time0
ensemble      = a->ensemble0

z!0 = "ensemble"
z!1 = "forecast_time"
z!2 = "lat"
z!3 = "lon"

z&lat           = lat
z&lon           = lon
z&forecast_time = forecast_time
z&ensemble      = ensemble

system("/bin/rm "+outfile)
fout = addfile(outfile,"c")
fout->Z0C  = z

end

 
 