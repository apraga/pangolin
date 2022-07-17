;*************************************************
; grid_fill_1.ncl
;*************************************************
;
; Concepts illustrated:
;   - Generating dummy data using "generate_2d_array"
;   - Using array syntax to arbitrarily set points to _FillValue
;   - Setting arguments for "poisson_grid_fill"
;   - Illustrating before and after application of poisson_grid_fill
;
;*************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
 ; nb_lat = 181 ;Contains both poles
  nb_lat = 180 
  nb_lat2 = 90 
  dlat = 1
  nb_lon = 1200 ;Less approximation error with a step of 0.1°
  ;nb_lon = 2*90+1
  dlon = 120./nb_lon
  nb_lon2 = 3*nb_lon

;******************************************************************
; Read file (contains only the concentration)
;******************************************************************
  filedir = "../plot/data/"
  filename = filedir + nb+".dat"

  nb_points = numAsciiRow(filename)
  data = asciiread(filename,(/nb_points,3/),"float")
  q = data(:,2)
  delete(data)

;******************************************************************
; Create 2D arrays for visu
;******************************************************************

  z  = new((/nb_lat,nb_lon2/),float)
  lat  = new((/nb_lat,nb_lon2/),float)
  lon  = new((/nb_lat,nb_lon2/),float)

  do i=0,nb_lat-1
    lat(i,0:nb_lon2-1) = tofloat(i+1)
  end do
  lat = 90. - lat*dlat

  do j=0,nb_lon2-1
    lon(0:nb_lat-1,j) = tofloat(j+1)
  end do
  lon = lon*dlon

  k = 0
  do i=1,nb_lat
    nb_mesh = 2*i-1
    if (i*dlat .gt. nb_lat2) then
      i2 = nb_lat+1 - i
      nb_mesh = 2*i2-1
    end if
    step_tmp = int2flt(nb_lon)/nb_mesh
    step = floattoint(step_tmp)
    nb_to_correct = nb_lon - nb_mesh*step
    dlon2 = 120./nb_mesh

    do nb_sec=0,2
      cur = nb_sec*nb_lon
      ;; No correction
      ;do j=0,nb_mesh-1
      ;  z(i,cur:cur+step-1) = rand()
      ;  cur = cur + step
      ;end do

      ; Correction of some cells sizes (we want to partition exactly a sector)
      do j=0,nb_to_correct-1
        z(i-1,cur:cur+step) = q(k)
        cur = cur + step+1
        k = k + 1
      end do
      do j=0,nb_mesh-nb_to_correct-1
        z(i-1,cur:cur+step-1) = q(k)
        cur = cur + step
        k = k + 1
      end do
    end do
    ;k = k + 2*nb_mesh
  end do

  ;wks   = gsn_open_wks("x11","raster")  
  output = "visu_raster_"+nb
   wks   = gsn_open_wks("ps",output)  
  gsn_define_colormap(wks,"BlueRed") ; "amwg")

  res = True
  ; Specify the mesh parameters
  res@sfXArray = lon
  res@sfYArray = lat

  res@cnFillOn = True            ; color Fill 
  res@cnFillMode = "CellFill"       ; faster than rasterfill
  res@gsnSpreadColors      = True                  ; use full range of colors
  res@cnLinesOn =  False             ; Turn off contour lines
  res@gsnMaximize = True    
  res@gsnAddCyclic = False
  res@cnMaxLevelValF  = 100.      ; set the maximum contour level

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

  plot = gsn_csm_contour_map(wks,z,res) 
end
