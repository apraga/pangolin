; Plot data with a given map projection for hdf5 format
; Data is 1D: q, lat, lon are each an 1D array

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
  ;folder = "/wkdir/pae2/praga/parallel/movie/gaussianhills_80lat_sp_CFL0.96/"
  folder = "../../generate/"
  date = "201301010000"
  file1 = addfile(folder + "ratio_1_" + date + ".h5","r")

  q = file1->ratio
 
  ;wks = gsn_open_wks("png", "files/"+date) ; Open a workstation.
  wks = gsn_open_wks("x11", "test")
  ;wks_type@wkOrientation = "landscape"
  res                     = True
  res@tiMainFuncCode    = ":"                 ; Needed for math mode
  ; Add a degree as a superscript
  res@tiMainString    = "Gaussian hills (1.13x0.67:F34:0:F:) "+date

  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@gsnSpreadColorStart = 24           ; at color index 24.

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode         = "AreaFill"   ; For contour
  ;res@cnFillMode         = "RasterFill"   ; For contour
                                       
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@sfXArray =  file1->lon
  res@sfYArray =  file1->lat

  ; Set levels manually
  res@cnLevelSelectionMode  = "ManualLevels"      ; Set the contour levels.
  res@cnMinLevelValF        =  0.0
  res@cnMaxLevelValF        = 1.1
  res@cnLevelSpacingF       =  0.01

  ; Using default projection cylindrical
  res@mpCenterLonF      =  180

  ; Polar
  ;res@gsnPolar = "NH"

  ;res@mpProjection = "WinkelTripel"
  ;res@mpProjection = "Mollweide"
;
;  map = gsn_csm_contour_map(wks,q(lev|0,lat|:,lon|:),res)
  map = gsn_csm_contour_map(wks,q,res)

end
