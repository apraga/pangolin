load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin

;******************************************************************
; Read file 
;******************************************************************
  filename = "data_orthographic_29112011.dat"
  nb_points = numAsciiRow(filename)
  data = asciiread(filename,(/nb_points,3/),"float")
  lat = data(:,0)
  lon = data(:,1)
  q = data(:,2)

  output = "visu_contour"
  wks   = gsn_open_wks("ps",output)  

  gsn_define_colormap(wks,"BlueRed")

  res = True
  ; Specify the mesh parameters
  res@sfXArray = lon
  res@sfYArray = lat

;*****************************************************************
; Color configuration
;*****************************************************************
  res@cnLinesOn=False
  res@cnFillOn = True 
  res@gsnSpreadColors = True ; use full range of colors

  res@gsnMaximize      = True    

  res@mpProjection = "Orthographic"
  res@mpCenterLatF = -1             ; Does not work
  res@mpCenterLonF = 20
  plot = gsn_csm_contour_map(wks,q,res) 

end

