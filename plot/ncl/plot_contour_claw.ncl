; Plot data with a given map pcv_raterojection

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
  ;folder = "/wkdir/pae2/praga/tmp/"
  ;folder = "/wkdir/pae2/praga/parallel/min_res/"
  ;name = "interp_0.75_0.out"
  name = "diff_claw_0.75.dat"
  ;name = str_concat((/folder, "ratio_1_201301122317.dat"/))
  ;name = str_concat((/folder, "ratio_800lat_un_CFL1.0_0.dat"/))

  ; Read file with unknown number of lines and split it
  print(name)
  data = asciiread(name, -1,"float")
  nrow = numAsciiRow(name)

  x = onedtond(data, (/nrow, 3/))
  pi = 4.*atan(1)
  q = x(:,2)
  lat = x(:,1)*180./pi
  lon = x(:,0)*180./pi

  ;wks = gsn_open_wks("x11","ratio_160lat_T")
  ;filen = str_concat((/"files/", date, ".png"/))
  wks = gsn_open_wks("png", "diff_claw_0.75.png")
  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnPaperOrientation = "portrait"
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnSpreadColorStart = 24           ; at color index 24.

  res@sfXArray            = lon
  res@sfYArray            = lat

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode         = "AreaFill"   ; For contour
                                       
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@cnLevelSelectionMode  = "ManualLevels"      ; Set the contour levels.
  res@cnMinLevelValF        =  0.0
  res@cnMaxLevelValF        = 1.0
  res@cnLevelSpacingF       =  0.02

  res@tiMainFuncCode    = ":"                 ; Needed for math mode
  ;stitle = str_concat((/date, " (0.56x0.37:F34:0:F:, 6 procs)"/))
  stitle = "Simu"
  res@tiMainString        = stitle

  ; Using default projection : cylindrical
  res@mpCenterLonF      =  0

  res@mpPerimOn         = False
  ; Disable continents fill 
  res@mpLandFillColor         = -1
  ; Disable outline
  ;res@mpOutlineOn = False;

  map = gsn_csm_contour_map(wks,q,res)
end
