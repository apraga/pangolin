; Plot data with a given map projection

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
  name = "/wkdir/pae2/praga/hourdin/ratio_160lat_CFL_0.7_T.dat"

  ; Read file with unknown number of lines and split it
  data = asciiread(name, -1,"float")
  nrow = numAsciiRow(name)

  x = onedtond(data, (/nrow, 3/))
  q = x(:,0)
  lat = x(:,1)
  lon = x(:,2)

  newlat = fspan(-90.,90,20)
  newlon = fspan(0.,360.,40)
  newlat@units = "degrees_north"
  newlon@units = "degrees_east"

  newq = dsgrid2(lon,lat,q,newlon,newlat)
  newlat@units = "degrees_north"
  newlon@units = "degrees_east"
  newq@lat = newlat
  newq@lon = newlon

  wks = gsn_open_wks("x11","ratio_160lat_T")
;  ;wks = gsn_open_wks("pdf","ratio_160lat_T")
  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnPaperOrientation = "portrait"
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnSpreadColorStart = 24           ; at color index 24.

;  res@sfXArray            = newlon
;  res@sfYArray            = newlat

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode         = "AreaFill"   ; For contour
                                       
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@cnLevelSelectionMode  = "ManualLevels"      ; Set the contour levels.
  res@cnMinLevelValF        =  0.0
  res@cnMaxLevelValF        = 1.0
  ;res@cnLevelSpacingF       =  0.02

  res@tiMainFuncCode    = ":"                 ; Needed for math mode
  res@tiMainString        = "Gaussian hills, T   (0.56x0.37:F34:0)"

  ; Using default projection : cylindrical
  res@mpCenterLonF      =  0

  res@mpPerimOn         = False
  ; Disable continents fill and outline
  res@mpLandFillColor         = -1
  res@mpOutlineOn = False;

  map = gsn_csm_contour_map(wks,newq,res)
end
