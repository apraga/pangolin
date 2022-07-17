
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
  nb_points = 100
  lat = new(nb_points,float)
  lon = new(nb_points,float)
  q = new(nb_points,float)
  dlat = 360./nb_points
  do i=0,nb_points/2-1
    lat(i) = 30.
    lat(2*i) = 20.
    lon(i) = 2*i*dlat 
    lon(2*i) = 2*i*dlat 
    q(i) = i
    q(2*i) = i
  end do

  wks   = gsn_open_wks("ps","projection")  

  res = True
  res@sfXArray = lon
  res@sfYArray = lat
  res@cnLinesOn=False
  res@cnFillOn = True 
  res@gsnMaximize = True    

  res@mpProjection = "Orthographic"
  ;res@mpCenterLatF = 30
  ;res@mpCenterLonF = 20

  plot = gsn_csm_contour_map(wks,q,res) 
end

