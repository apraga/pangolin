; Set the map resources as Lauritzen

  undef("setRes")
procedure setResLauritzen (res: logical)
  begin
  res@mpGridAndLimbOn = True
  res@mpGridSpacingF = 90.0
  res@mpGridLineDashPattern = 2
  res@mpLimbLineColor = 1
  res@mpGridLineColor = 0
  res@mpGreatCircleLinesOn = True
  res@mpPerimLineThicknessF   = 1.

  res@mpCenterLonF =  180.0           ; Rotate the projection.
  res@mpOutlineOn     = False       ; No continents

  mnlvl = 0.0                 ; Minimum contour level.
  mxlvl = 1.1               ; Maximum contour level.
  spcng = 0.05                  ; Contour level spacing.
  ncn   = (mxlvl-mnlvl)/spcng + 1  ; Number of contour levels.

  ;
  ; Contour line features
  ;

  res@cnLevelSelectionMode = "ManualLevels" ; Define your own
  res@cnMinLevelValF       = mnlvl          ; contour levels.
  res@cnMaxLevelValF       = mxlvl
  res@cnLevelSpacingF      = spcng

  res@cnInfoLabelOn        = False    ;No c-label box at the bottom 
  res@cnLineLabelsOn       = False    ;No cotour labels 
  res@cnInfoLabelPerimOn   = False

  res@cnFillOn        = True               ; color plot desired
  res@gsnSpreadColors = True               ; use full color table
  res@gsnSpreadColorStart = 10
  res@gsnSpreadColorEnd   = -4

  res@lbLabelBarOn     = False

  res@lbLabelStride       = 2       ;
  res@lbTitleFontThicknessF   = .08
  res@lbAutoManage = False
  res@lbLabelFontHeightF    = .0175

  ; Default ticks make weird plot with columns and slotted cylinders
  ;res@pmTickMarkDisplayMode = "Always"
  ; We redefine them
  res@tmXBMinorOn = False
  res@tmYLMinorOn = False 
  res@gsnMajorLonSpacing = 90   
  res@gsnMajorLatSpacing = 45   

  res@cnLinesOn         = False
  res@cnLineDashPattern    = 0     ; Solid line contours     
  res@cnLineThicknessF     = 1.0   ; Double the line thickness.
  res@cnMonoLineColor      = True
  res@cnLineColors         = 1 

  res@cnLevelFlag       = True
  res@cnFillOn          = True

  res@mpPerimOn         = False
  res@mpLandFillColor         = -1 ; Disable continents fill 
  ; Disable outline
  res@mpOutlineOn = False;

  ;fontheight = .027
  ;res@txFont = "helvetica"
  ;res@gsnStringFontHeightF = fontheight

  end

; Set map resources for columns plot
procedure setPanelResLauritzen (pres: logical, titlelist: string)
  begin

  whiteSpace = " "
  pres@txString = "         "+"0.75~F34~0~F~, unlimited            "+\
                   whiteSpace+"  "+ "0.75~F34~0~F~, shape-preserving"

  pres@txFontHeightF = 0.013

  pres@gsnPanelBottom   = 0.05
  pres@gsnPanelLabelBar  = True
  pres@lbLabelAutoStride = True

  pres@gsnPanelXWhiteSpacePercent =5 
  pres@gsnPanelYWhiteSpacePercent =7

  pres@amJust = "BottomLeft"
  pres@gsnPanelFigureStrings= titlelist
  pres@gsnPanelFigureStringsFontHeightF = 0.008
  pres@gsnPanelLabelBar  = True
  pres@lbLabelFontHeightF = 0.012
  pres@lbLabelStride = 2

  end
