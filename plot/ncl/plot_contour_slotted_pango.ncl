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
if (j .eq. 0) then
  if (i .eq. 0) then
    title = "Unlimited"
  else
    title = "Shape-preserving"
  end if
  ;title = str_concat((/title, " (0.56x0.37:F34:0:F:)"/))
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
  res@tiYAxisString = "t=" + time2
end if


end


begin

wks = gsn_open_wks("pdf", "slotted_cylinders_pango.pdf")
res                     = True
res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
gsn_define_colormap(wks,"wh-bl-gr-ye-re")

; Lauritzen config
setResLauritzen(res)


times = (/"0", "T"/)
types = (/"un", "sp"/)
folder = "/wkdir/pae2/praga/paper/"
rows = dimsizes(times)
cols = dimsizes(types)

map = new(rows*cols, graphic) 
k = 0
do j = 0, rows-1
  do i = 0, cols-1
    name = folder + "ratio_160lat_slotted_" + types(i) + "_CFL0.96_" + times(j) + ".dat"
    set_title(res, i, j, times(j))
    pangol_ascii_contour(wks, res, map, name, k)
    k = k + 1
  end do
end do

;************************************************
; create panel
;************************************************
resP                     = True                ; modify the panel plot
;resP@gsnPanelXWhiteSpacePercent =50
;resP@gsnPanelYWhiteSpacePercent =7

resP@tiMainFuncCode    = ":"                 ; Needed for math mode
title = "Slotted cylinders, t=0, t=T"
title = str_concat((/title, " (0.56x0.37:F34:0:F:)"/))
;resP@txString            = title
resP@gsnPanelLabelBar    = True                ; add common colorbar
gsn_panel(wks, map, (/rows,cols/), resP)               ; now draw as one plot

end
