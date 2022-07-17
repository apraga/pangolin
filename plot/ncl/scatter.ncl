; Plot data with a given map projection

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
  name = "diff.dat";

  ; Read file with unknown number of lines and split it
  data = asciiread(name, -1,"float")
  nrow = numAsciiRow(name)

  x = onedtond(data, (/nrow, 3/))
  q = x(:,0)
  lat = x(:,1)
  lon = x(:,2)

  wks = gsn_open_wks("x11","ratio_160lat_T")
  ;filen = str_concat((/"files/", date, ".png"/))
  ;wks = gsn_open_wks("png", "un.png")
  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
                                       
  res@gsnDraw            = False        ; Don't draw plot or advance
  res@gsnFrame           = False        ; frame. Will do this later.

  res@tiMainString       = "Scatter plot with grouped markers"

  levels  = ispan(0,10,1)*0.1
  nlevels = dimsizes(levels)

  nstep = 256/nlevels
  colors = ispan(2,256,nstep)
;---Set some marker resources
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerThicknessF = 2.5
  res@xyMarkerColors     = colors
  ;res@xyMarker           = square     ; this is a filled square

;  labels  = new(nlevels,string)
;  labels  = ""
;  labels(0::10) = "" + ispan(0,50,10)
  
  print("min/max data = " + min(q) + "/" + max(q))
  print("# of groups  = " + (nlevels-1))
  npts = dimsizes(q)

;---Create new 2D array to hold groupings of values
  data2d = new((/nlevels-1,npts/),typeof(q))

;---Group the values and put in 2D array.
  do i=0,nlevels-2
    ii = ind(q.ge.levels(i).and.q.lt.levels(i+1))
    ;Cif (.not. ismissing(ii)) then
    if (.not. all(ismissing(ii))) then
      data2d(i,ii) = q(ii)
    end if
    delete(ii)
  end do
  plot  = gsn_csm_y (wks,data2d,res)
  maximize_output(wks,False)
end
