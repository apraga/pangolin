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
  nb_lat = 181 ;Contains both poles
  dlat = 1
  nb_lon = 1200 ;Less approximation error with a step of 0.1°
  ;nb_lon = 2*90+1
  dlon = 120./nb_lon
  nb_lon2 = 3*nb_lon

;******************************************************************
; Read file (contains only the concentration)
;******************************************************************
  filedir = "../plot/data/"
  filename = filedir + "0"+".dat"

  nb_points = numAsciiRow(filename)
  data = asciiread(filename,nb_points,"float")

;******************************************************************
; Create 2D arrays for visu
;******************************************************************

  z  = new((/nb_lat,nb_lon2/),float)
  lat  = new((/nb_lat,nb_lon2/),float)
  lon  = new((/nb_lat,nb_lon2/),float)

  do i=0,nb_lat-1
    lat(i,0:nb_lon2-1) = tofloat(i)
  end do
  lat = 90. - lat*dlat

  do j=0,nb_lon2-1
    lon(0:nb_lat-1,j) = tofloat(j)
  end do
  lon = lon*dlon

  k = 0
  do i=0,nb_lat-1
    nb_mesh = 2*i+1
    if (i .gt. 90) then
      nb_mesh = 2*(nb_lat-i-1)-1
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
        z(i,cur:cur+step) = data(k);rand()
        cur = cur + step+1
        k = k + 1
      end do
      do j=0,nb_mesh-nb_to_correct-1
        ;print(sprintf("%f",k)+sprintf("%f",i))
        z(i,cur:cur+step-1) = data(k);and()
        cur = cur + step
        k = k + 1
      end do
    end do
  end do

  wks   = gsn_open_wks("x11","raster")  
  ; wks   = gsn_open_wks("ps","raster_demo")  
  ;gsn_define_colormap(wks,"gui_default") ; "amwg")

  res = True
  ; Specify the mesh parameters
  res@sfXArray = lon
  res@sfYArray = lat

  res@cnFillOn = True            ; color Fill 
  res@cnFillMode = "CellFill"       ; faster than rasterfill
  res@cnLinesOn =  False             ; Turn off contour lines
  res@gsnMaximize = True    
  res@gsnAddCyclic = False

  res@mpProjection          = "Orthographic"
  res@mpCenterLatF          =  20
  res@mpCenterLonF          =  90

  plot = gsn_csm_contour_map(wks,z,res) 
end
