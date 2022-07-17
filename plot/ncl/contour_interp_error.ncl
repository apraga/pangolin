; Plot data with a given map projection

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
  title = "4x4:F34:0:F: lat-lon to 2x2:F34:0:F: lat-lon, bilinear"
  output = "4x4_reg_2x2_reg_conservative"
  projection = "WinkelTripel"
  ;plot = "normal";error"
  plot = "error"


  folder = "../../netcdf/"
  file2 = addfile(folder + "fine_interp.nc","r")

  if (plot .eq. "error") then
    file1 = addfile(folder + "fine_init.nc","r")
    q = file2->temperature - file1->temperature
  else
    q = file2->temperature
  end if
 
  q!1    = "lat"          
  q!2    = "lon"        
  q!0    = "lev"        
  q&lat  = file2->lat    
  q&lon  = file2->lon  
  q&lev  = file2->lev  
  q@nlon = dimsizes(q&lon)
  q@nlat = dimsizes(q&lat) 
 
  wks_type = "pdf"
  ;wks_type = "x11"
  wks = gsn_open_wks(wks_type, output)
  res                     = True
  res@tiMainFuncCode    = ":"                 ; Needed for math mode
  res@tiMainString    = title

  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnSpreadColorStart = 24           ; at color index 24.

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode         = "AreaFill"   ; For contour
                                       
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  ; Log plot from the max and min
  if (plot .eq. "error") then
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    levls = new(15, double)
    n = dimsizes(levls)
    levls(0) = min(q)
    levls(n-1) = max(q)
    do i = 1, n/2-1
      levls(i) = levls(i-1)*0.1
      levls(n-1-i) = levls(n-1-(i-1))*0.1
    end do
    res@cnLevels = levls ; set levels 
  else
    res@cnLevelSelectionMode  = "ManualLevels"      ; Set the contour levels.
    res@cnMinLevelValF        =  0.0
    res@cnMaxLevelValF        = 1.0
    res@cnLevelSpacingF       =  0.02
  end if


  if (projection .eq. "polar") then
    res@gsnPolar = "SH"
  else
    res@mpProjection = projection
  end if

;
  map = gsn_csm_contour_map(wks,q(lev|0,lat|:,lon|:),res)

end
