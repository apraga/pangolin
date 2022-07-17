; Read data, create contour and put it into and array
procedure pangol_ascii_contour(wks: graphic, res: logical, map: graphic, name: string, k: integer)
begin
  ; Read file with unknown number of lines and split it
  ;print(name)
  data = asciiread(name, -1,"float")
  nrow = numAsciiRow(name)

  x = onedtond(data, (/nrow, 3/))
  q = x(:,0)
  lat = x(:,1)
  lon = x(:,2)

  res@sfXArray            = lon
  res@sfYArray            = lat
  map(k) = gsn_csm_contour_map(wks,q,res)

  delete(res@sfXArray)      
  delete(res@sfYArray)      
  delete(q)
  delete(lat)
  delete(lon)

end

; Store plot in graphic array
; f is given by the addfile function
procedure lauritzen_contour(wks: graphic, res: logical, map: graphic, \
    name: string, k: integer)
begin
  ;print(name)
  f = addfile(name, "r")
  n = getfilevardimsizes(f, "data") 
  lon = f->lon
  lat = f->lat

  if (dimsizes(n) .eq. 3) then
    data = f ->data(0,:,:)

  else
    if (dimsizes(n) .eq. 2) then
      data = f ->data(:,:)
      ; Attach coordinates data
      data!0 = "lat"
      data!1 = "lon"
      data&lat = lat
      data&lon = lon
      res@gsnAddCyclic = False                       ; Otherwise, wrong coordinates
    else
      data = f->data(:)
      res@sfXArray   = lon     
      res@sfYArray   = lat    

    end if
  end if

  res@mpCenterLonF = 180. 

  res@gsnLeftString   = ""
  res@gsnCenterString = ""
  res@gsnRightString  = ""

  map(k) = gsn_csm_contour_map_ce(wks,data,res) 
  if (dimsizes(n) .eq. 1) then
    delete(res@sfXArray)      
    delete(res@sfYArray)      
  end if
  delete(res@mpCenterLonF)

  delete(data)
  delete(lat)
  delete(lon)

  end
