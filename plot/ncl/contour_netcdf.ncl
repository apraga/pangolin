; Plot data with a given map projection

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
  folder = "../../netcdf/"
  file1 = addfile(folder + "fine_init.nc","r")
  file2 = addfile(folder + "fine_interp.nc","r")

  q1 = file1->temperature
  q2 = file2->temperature
  q = q2 - q1
 
  q!1    = "lat"          
  q!2    = "lon"        
  q!0    = "lev"        
  q&lat  = file1->lat    
  q&lon  = file1->lon  
  q&lev  = file1->lev  
  q@nlon = dimsizes(q&lon)
  q@nlat = dimsizes(q&lat) 
 
  ;!wks = gsn_open_wks("x11","gsun05n") ; Open a workstation.
  wks_type = "png"
  ;wks_type = "x11"
  ;wks_type@wkOrientation = "landscape"
  wks = gsn_open_wks(wks_type, "fine_diff");fine_to_coarse")
  res                     = True
  res@tiMainFuncCode    = ":"                 ; Needed for math mode
  res@tiMainString    = "Conservative interpolation (2x2:F34:0:F: lat-lon to 4x4:F34:0:F: lat-lon )"

  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnSpreadColorStart = 24           ; at color index 24.

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode         = "AreaFill"   ; For contour
  ;res@cnFillMode         = "RasterFill"   ; For contour
                                       
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  ; Set the number of levels
  ;res@cnLevelSelectionMode  = "EqualSpacedLevels"      ; Set the contour levels.
  ;res@cnMaxLevelCount  = 100

  ; Log plot from the max and min
  res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
  levls = new(15, double)
  n = dimsizes(levls)
  levls(0) = -0.03;min(q)
  levls(n-1) = 2;max(q)
  do i = 1, n/2-1
    levls(i) = levls(i-1)*0.1
    levls(n-1-i) = levls(n-1-(i-1))*0.1
  end do
  print(levls)
  res@cnLevels = levls ; set levels 

  ;; Set levels manuall
  ;res@cnLevelSelectionMode  = "ManualLevels"      ; Set the contour levels.
  ;res@cnMinLevelValF        =  -1.0
  ;res@cnMaxLevelValF        = 1.0
  ;res@cnLevelSpacingF       =  0.1

  ; Using default projection : cylindrical
  res@mpCenterLonF      =  180

  ; Polar
  ;res@gsnPolar = "NH"

  ;res@mpProjection = "WinkelTripel"
  res@mpProjection = "Mollweide"
;
  map = gsn_csm_contour_map(wks,q(lev|0,lat|:,lon|:),res)

end
