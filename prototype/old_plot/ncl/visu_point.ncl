;********************************************************
; Create scattered plot for irregular grid              ;
;********************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  ;********************************************************
  ;                 Read a file                        ;
  ;********************************************************
  filedir = "../plot/data/"
  ;filedir = "../plot/data_hourdin/"
  filename = filedir + nb+".dat"

  nb_points = numAsciiRow(filename)
  q = asciiread(filename,(/nb_points,3/),"float")
  markers = new((/nb_points/),graphic)

  ;********************************************************
  ;              Plotting parameters                   ;
  ;********************************************************
  ;wks   = gsn_open_wks("x11",type)  
  wks   = gsn_open_wks("ps","visu_point_"+nb)  
  res = True
  ; Type of projection
  if (type .eq. "ortho") then ; Center of projection
    res@mpProjection = "Orthographic"
    res@mpCenterLatF = 0
    ;res@mpCenterLonF = 90
  else
    if (type .eq. "lambert") then ; Center of projection
      res@mpProjection = "LambertConformal"
      res@mpLimitMode = "Corners"  
      res@mpLeftCornerLatF = 60
      res@mpLeftCornerLonF = -30
      res@mpRightCornerLatF = 30 
      res@mpRightCornerLonF = 30 
    else
      if (type .eq. "satellite") then ; Center of projection
        res@mpProjection = "Satellite" ; choose map projection
        res@mpCenterLonF = 0.0       ; choose center lon
        res@mpCenterLatF = 40.         ; choose center lat
        res@mpSatelliteDistF = 5.0         ; choose satellite
      else
  ;      res@mpCenterLonF = 90.  
      end if
    end if
  end if

  res@gsnFrame    = False                             ; don't advance yet
  res@gsnMaximize = True
  ; Get color table
  gsn_define_colormap(wks,"rainbow")
  rgb_array = gsn_retrieve_colormap(wks)
  size_colors = dimsizes(rgb_array)
  nb_colors = size_colors(0)

  ; Special colors for oceans and lands
  dark = NhlNewColor(wks,0.,0.,0.)
  grey_light = NhlNewColor(wks,0.8,0.8,0.8)
  res@mpOceanFillColor = grey_light 
  res@mpLandFillColor = dark
  res@mpInlandWaterFillColor = grey_light
 
  ; If ortho, place it by hand
  if (type .eq. "ortho") then 
    res@gsnMaximize = False
    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.8
    res@vpHeightF = 0.8
    res@vpXF= 0.18
    res@vpYF= 0.9
  end if
  map = gsn_csm_map(wks,res)                ; Create the map plot, don't draw
  
  ;****************************************************
  ;              Create markers                        ;
  ;****************************************************
  mkres = True
  ; Size of the markers
  mkres@gsMarkerIndex = 16                
  ; mkres@gsMarkerSizeF = 0.01 ; Without holes
  ;mkres@gsMarkerSizeF = 0.0025 ; Too small
  mkres@gsMarkerSizeF = 0.005 ; Good compromise

  mkres@pmLabelBarDisplayMode = "Always"          ;

  c = q(:,2)
  ;q_max = 181.;max(c)
  q_max = max(c)
  q_min = min(c)

  step_lat = 1 ; Skip some latitudes if > 1

  ; Regroup data by values (faster)
  levels  = fspan(q_min,q_max,nb_colors-1)

  ;---Create new 2D array to hold groupings of values 
  lat = new((/nb_colors-1,nb_points/),float)
  lon = new((/nb_colors-1,nb_points/),float)

  do i=2,nb_colors-3 ;Faster without the 0
  ;do i=0,nb_colors-3
    ii = ind(c.ge.levels(i).and.c.lt.levels(i+1))
    size = dimsizes(ii)
    if (size(0) .lt. 2) then
      delete(ii)
      continue
    end if
    lat(i,ii) = q(ii,0)
    lon(i,ii) = q(ii,1)
    delete(ii)
    mkres@gsMarkerColor = i+2 ; No white
    gsn_polymarker(wks,map,lon(i,:),lat(i,:),mkres)
  end do

  ; Add the missing values
  ii = ind(c.ge.levels(nb_colors-2))
  lat(nb_colors-2,ii) = q(ii,0)
  lon(nb_colors-2,ii) = q(ii,1)
  delete(ii)
  mkres@gsMarkerColor = nb_colors-2 ; No white
  gsn_polymarker(wks,map,lon(nb_colors-2,:),lat(nb_colors-2,:),mkres)


  ; ******************************************************** ;
  ;                         Label bar
  ; ******************************************************** ;
  delta_q = q_max-q_min
  if (delta_q .le. 1.0) then
    ; Skip white colors
    step = delta_q/nb_colors
    start_q = q_min+2*step
    tmp = fspan(start_q,q_max,nb_colors-2)
    labels = " " + tofloat(tmp)
  else
    tmp = ispan(2,nb_colors-1,1)*delta_q/nb_colors   ; Labels for boxes.
    tmp = q_min + tmp
    labels = " " + toint(tmp)
  end if
  
  lbres = True
  lbres@lbAutoManage = False   ; Allows us to set more 
  lbres@lbTitleFontHeightF= 0.02    ; Font heights for main title
  lbres@lbLabelFontHeightF= 0.015   ; and box labels

  if (type .eq. "ortho") then ; Center of projection
    lbres@vpWidthF = 0.20    ; Width of labelbar
    lbres@vpHeightF = 1.0    ; Height of labelbar
    pos_x = 0.0 
    pos_y = 1.0 
  else
    lbres@lbOrientation = "Horizontal"
    lbres@lbTitleFontHeightF= 0.02    ; Font heights for main title
    lbres@lbLabelFontHeightF= 0.015   ; and box labels
    lbres@vpWidthF = 0.80    ; Width of labelbar
    lbres@vpHeightF = 0.20    ; Height of labelbar
    pos_x = 0.1 
    pos_y = 0.27 
  end if

  lbres@lbLabelAutoStride = True    ; Controls label spacing so labels
  lbres@lbMonoFillPattern = True
  lbres@lbFillColors = ispan(2,nb_colors-1,1)
  lbres@lbPerimOn = False   ; Turn off box around labelbar.
  lbres@lbBoxLinesOn= "False"

  lbres@lbTitleString = "Concentration"
  gsn_labelbar_ndc(wks,nb_colors,labels,pos_x,pos_y,lbres)  ; Draw a labelbar

  frame(wks)
end
