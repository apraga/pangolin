load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************************
begin
  name = "../../data/filaments/filaments_limited/ratio_80lat_T2.dat"
  ; Read file with unknown number of lines and split it
  data = asciiread(name, -1,"float")
  nrow = numAsciiRow(name)

  x = onedtond(data, (/nrow, 3/))
  q = x(:,0)
  lat = x(:,1)
  lon = x(:,2)

  wks = gsn_open_wks("x11","lol")
  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnSpreadColorStart = 24           ; at color index 24.

  res@sfXArray            = lon
  res@sfYArray            = lat

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Only one correct
                                       
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@cnLevelSelectionMode  = "ManualLevels"      ; Set the contour levels.
  res@cnMinLevelValF        =  0.0
  res@cnMaxLevelValF        = 1.0
  res@cnLevelSpacingF       =  0.01

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.

  res@tiMainString        = "Pangolin"


  ;res@mpProjection      = "Mollweide" ; Projection
  ;res@mpCenterLatF      = 0
  ;res@mpCenterLonF      =  0

  res@mpPerimOn         = False

  map = gsn_csm_contour_map(wks,q,res)     ; Draw a contour plot.
end

