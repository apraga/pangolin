begin

; Initialize variables
projection = (/"orthographic","mollweide","stereographic"/)

;
; Create Application object.  The Application object name is used to
; determine the name of the resource file, which is "pr05.res" in this
; case.
; 
appid = create "pr05" appClass defaultapp 
    "appDefaultParent" : True
    "appUsrDir" : "./"
end create

;
; Default is to display output to an X workstation
;
X11=0
PS=1

if (X11 .eq. 1) then
  wid = create "pr05Work" xWorkstationClass defaultapp
      "wkPause" : True
  end create
else
  if (PS .eq. 1) then
    output = "visu_polygon_oo_" + nb+".ps"
    wid = create "pr05Work" psWorkstationClass defaultapp
      "wkColorMap" : "BlueRed"
      "wkPSFileName" : output
    end create
  end if
end if

;
; Create a MapPlot that covers the entire NDC space 
; to use as a drawing canvas
;
canvas = create "canvas" mapPlotClass wid
  "vpXF" : -0.3
  "vpYF" : 1.3
  "vpWidthF" : 1.5
  "vpHeightF" : 1.5
end create

;
; Create a GraphicStyle to control the primitive attributes.
;
gsid = create "style" graphicStyleClass wid end create

i = 0
;
; Set the map projection
;
  setvalues canvas
    "mpProjection" : projection(0)
    "mpCenterLatF" : 30
    "mpCenterLonF" : 50

  end setvalues
;
  setvalues gsid
    "gsMarkerColor" : 0
    "gsEdgesOn" : False
  end setvalues

  plon = new (4,float)
  plat = new (4,float)

  nb_lat = 179 
  nb_lat2 = 90 
  dlat = 1
  nb_lon = 2*90+1
  dlon = 120./nb_lon
  nb_lon2 = 3*nb_lon

;******************************************************************
; Read file (contains only the concentration)
;******************************************************************
  filedir = "../plot/data/"
  filename = filedir + nb+".dat"

  nb_points = 3*nb_lat2*nb_lat2 + 3*(nb_lat2-1)*(nb_lat2-1); numAsciiRow(filename)
  data = asciiread(filename,(/nb_points,3/),"float")
  q = data(:,2)

  dlat = 1
  j = 0
  k = 0

  q_max = max(q)
  q_min = min(q)
  delta_q = q_max - q_min

  do i=0,nb_lat-1
    nb_mesh = 3*(2*i - 1)
    if (i .gt. 90) then
      i2 = nb_lat+1-i
      nb_mesh = 3*(2*i2-1)
    end if
    dlon = 360./nb_mesh
    do j=0,nb_mesh-1
      value = (q(k) - q_min)/delta_q
      color = value*251+2
      setvalues gsid
        "gsFillColor" : color
      end setvalues
      k = k + 1

      lon = j*dlon
      lat = 89-i*dlat
      plon(0) = lon
      plon(1) = lon+dlon
      plon(2) = lon+dlon
      plon(3) = lon

      plat(0) = lat
      plat(1) = lat
      plat(2) = lat+dlat
      plat(3) = lat+dlat
      NhlDataPolygon(canvas,gsid,plon,plat)
    end do
  end do

; Draw the map outlines and grid
  draw(canvas)
  frame(wid)

end
