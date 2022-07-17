load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
  wks   = gsn_open_wks("x11","traj")  
  nb_lat = 90
  nb_lon = 360
  lat = new ((/nb_lat/),float)
  lon = new ((/nb_lon/),float)
  q = new ((/nb_lat,nb_lon/),float)
  dlat = 1.
  dlon = 1.

  do i=0,nb_lat-1
    lat(i) = i*dlat
  end do
  do j=0,nb_lon-1
    lon(j) = j*dlon
  end do

  do i=1,nb_lat
    do j=1,nb_lon
      cos_lat = cos(i*dlat)
      cos_lon = cos(j*dlon*0.5)
      q(i-1,j-1) = cos_lat*cos_lat*cos_lon*cos_lon*cos_lon
    end do
  end do
  q!0    = "lat"                 ; Name dimensions 0 and 1 of sst
  q!1    = "lon"                 ; "lat" and "lon.
  q&lat  = lat
  q&lon  = lon
  q&lat@units = "degrees_north"  ; Attach units attributes
  q&lon@units = "degrees_east"
  res             = True
  res@mpFillOn    = True
  res@cnLevelSelectionMode = "ManualLevels" ; manually set the contour levels
    res@cnMinLevelValF  = -.9     ; set the minimum contour level
      res@cnMaxLevelValF  = .9     ; set the maximum contour level
        res@cnLevelSpacingF = .5      ; set the interval between contours
  ;res@sfXArray = lat
  ;res@sfYArray = lon
  plot = gsn_csm_contour_map(wks,q,res)  

end

