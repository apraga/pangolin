load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin

;******************************************************************
; Read file (contains the concentration, and lat-lon)
;******************************************************************
  filedir = "../plot/data/"
  ;filedir = "../plot/data_old/"
  filename = filedir + nb+".dat"

  nb_points = numAsciiRow(filename)
  data = asciiread(filename,(/nb_points,3/),"float")
  lat = data(:,0)
  lon = data(:,1)
  q = data(:,2)
  print(max(q))
  print(min(q))

  ;wks   = gsn_open_wks("x11","raster")  
  output = "visu_contour_"+nb
  wks   = gsn_open_wks("ps",output)  

  ;gsn_define_colormap(wks,"gui_default") ; "amwg")
  gsn_define_colormap(wks,"BlueRed")

  res = True
  ; Specify the mesh parameters
  res@sfXArray = lon
  res@sfYArray = lat
;  res@trGridType          = "TriangularMesh"    ; Define the mesh type.  

;*****************************************************************
; Color configuration
;*****************************************************************
  res@cnLinesOn=False
  res@cnFillOn         = True            ; color Fill 
  ;res@cnFillMode       = "CellFill"       ; Cellfill mode

  res@gsnSpreadColors      = True                  ; use full range of colors

  res@gsnMaximize      = True    
  ;res@gsnAddCyclic=True

  ;res@cnLevelSelectionMode = "ManualLevels" ; manually set the contour levels
  ;res@cnMinLevelValF  = 0.     ; set the minimum contour level
  ;res@cnMaxLevelValF  = 100.     ; set the maximum contour level
  ;res@cnLevelSpacingF = 1.      ; set the interval between contours

  if (type .eq. "ortho") then ; Center of projection
    res@mpProjection = "Orthographic"
    res@mpCenterLatF = -80
    res@mpCenterLonF = 20
    plot = gsn_csm_contour_map(wks,q,res) 
  else
    if (type .eq. "polar") then ; Center of projection
      res@gsnPolar   = "SH"                          ; specify the hemisphere
      plot = gsn_csm_contour_map_polar(wks,q,res)    ; create the plot
    else
      res@mpCenterLonF = 90
      plot = gsn_csm_contour_map(wks,q,res) 
    end if
  end if


end

