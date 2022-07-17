
;*************************************************
; polyg_4.ncl
;
; Concepts illustrated:
;   - Drawing a cylindrical equidistant map
;   - Zooming in on a particular area on a cylindrical equidistant map
;   - Attaching an outlined box to a map plot
;   - Attaching filled polygons to a map plot
;   - Filling in polygons with a shaded pattern
;   - Changing the color and thickness of polylines
;   - Changing the color of a filled polygon
;   - Labeling the lines in a polyline
;   - Changing the density of a fill pattern
;   - Adding text to a plot
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
;************************************************
begin
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","polyg")                ; open a ncgm file

  res           = True                              ; plot mods desired
  res@gsnDraw   = False                             ; don't draw yet
  res@gsnFrame  = False                             ; don't advance yet

  ;plot = gsn_map(wks,"Orthographic",res) ; create the plot   
  ;plot = gsn_map(wks,"Mollweide",res) ; create the plot   
  plot = gsn_csm_map_ce(wks,res) ; create the plot   

  
  draw(plot)
;************************************************
; add the box
;************************************************
  resp                  = True                      ; polyline mods desired

; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.

  nb_lat2 = 90
  nb_lat = 2*nb_lat2
  xcur = new(4, float)
  ycur = new(4, float)
  j = 40
  dlat = 90./nb_lat2

  srand (12345)     
  do i = 1, nb_lat
    i2 = i
    if (i .gt. nb_lat2) then
      i2 = nb_lat - i + 1
    end if
    nb_lon = 3*(2*i2-1)
    dlon = 360./nb_lon
    lat = 90.-(i-1)*dlat
    do j = 0, nb_lon-1
      rand_num = (rand()/128) + 1      
      resp@gsFillColor = rand_num 
      xcur = (/j*dlon, j*dlon, (j+1)*dlon, (j+1)*dlon/)
      ycur = (/lat - dlat, lat, lat, lat -dlat/)
      gsn_polygon(wks,plot,xcur,ycur,resp)      
    end do
  end do


  frame(wks)

end
