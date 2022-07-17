% Write data on the reduced grid
% Format : (latitude, longitude, value)
% Latitude range : 0:360
% Longitude range : 90:-90 degree
% Cells are browsed at constant latitude
function write_data(k,x,nb_lat,nb_lat2,dlat)
  name = strcat('data/',num2str(k),'.dat');
  fid= fopen (name,'w');
  precision = '16'; % Precision machine
  format = strcat('%.',precision,'g %.',precision,'g %.',precision,'g \n');
  if (fid == -1)
    disp('error opening file ');
  end

  for i = 1:nb_lat
    nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
    nb_cells = 3*nb_cells;
    lat = 90-(i-0.5)*dlat;
    dlon = 360/(nb_cells);
    for j=1:nb_cells
      lon = (j-0.5)*dlon;
      count= fprintf (fid,format, x(i,j), lat, lon);
    end
  end
  fclose(fid);
end
