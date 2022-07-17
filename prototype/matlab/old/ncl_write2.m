% Ecrit dans un fichier des données sur une grille reguliere 
%
function ncl_write2(k,x,nb_lat,nb_lon)
  name = strcat('../../plot/data_old/',num2str(k),'.dat')
  fid= fopen (name,'w');
  if (fid == -1)
    disp('error opening file ')
  end

  dlat = 1;
  k = 1;

  for i = 1:nb_lat
    nb_mesh = 3*(2*i-1);
    lat = 90-i*dlat;
    dlon = 360/(nb_mesh);
    for j=1:nb_mesh
      lon = j*dlon;
      count= fprintf (fid,'%f %f %f \n',lat,lon,x(k));
      k = k + 1;
    end
  end
  fclose(fid);
end
