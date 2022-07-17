% Ecrit dans un fichier des données sur une grille reguliere 
% Latitude : 0:360
% Longitude : 90:-90 degree
function ncl_write(k,x,nb_lat,nb_lat2,nb_lon)
  name = strcat('../plot/data/',num2str(k),'.dat');
  fid= fopen (name,'w');
  if (fid == -1)
    disp('error opening file ');
  end

  dlat = 1;

  for i = 1:nb_lat
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    nb_mesh = 3*nb_mesh;
    lat = 90.5-i*dlat;
    dlon = 360/(nb_mesh);
    for j=1:nb_mesh
      lon = (j-0.5)*dlon;
      count= fprintf (fid,'%f %f %f \n',lat,lon,x(i,j));
    end
  end
  fclose(fid);
end
