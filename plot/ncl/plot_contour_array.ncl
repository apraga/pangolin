; Array of plots with a given map pcv_raterojection
; PLots are identical to the ones in Lauritzen2014

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "set_res_lauritzen.ncl"
load "one_contour_plot.ncl"

; Must be called before one_plot
procedure set_title(res: logical, i: integer, j: integer, time: string)
begin
; Main title for the column
if (j .eq. 0) then
  if (i .eq. 0) then
    title = "Cosine bells"
  else
    title = "Gaussian hills"
  end if
  title = str_concat((/title, " (0.56x0.37:F34:0:F:)"/))
  res@tiMainFuncCode    = ":"                 ; Needed for math mode
  res@tiMainString = title
else
  res@tiMainOn = False
end if

; Left axis title
time2 = time
res@tiYAxisOn = False
if (i .eq. 0) then
  res@tiYAxisOn = True
  res@tiYAxisFontHeightF = 0.04
  if (time2 .eq. "T2") then
    time2 = "T/2"
  end if
  res@tiYAxisString = "t=" + time2
end if

end

;************************************************
begin

wks = gsn_open_wks("pdf", "gaussian_cosine_results.pdf")
res                     = True
res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame

; Lauritzen config
setResLauritzen(res)

map = new(6, graphic)                          ; create a plot array

times = (/"0", "T2", "T"/)
types = (/"cosine", "gaussian"/)
folder = "/wkdir/pae2/praga/paper/"

k = 0
do j = 0, dimsizes(times)-1
  do i = 0, dimsizes(types)-1
    name = folder + "ratio_160lat_" + types(i) + "_sp_CFL0.96_" + times(j) + ".dat"
    set_title(res, i, j, times(j))
    pangol_ascii_contour(wks, res, map, name, k)

    k = k + 1
  end do
end do

;************************************************
; create panel
;************************************************
resP                     = True                ; modify the panel plot
resP@tiMainFuncCode    = ":"                 ; Needed for math mode
title = "Cosine bells and Gaussian hills, t=0, t=T/2, t=T"
title = str_concat((/title, " (0.56x0.37:F34:0:F:)"/))
;resP@txString            = title
resP@gsnPanelLabelBar    = True                ; add common colorbar
gsn_panel(wks, map, (/3,2/), resP)               ; now draw as one plot

end
