function ncl_write(u_s,v_s,nb_lat,nb_lat2,nb_lon)
  name = strcat('../../plot/data/data_vectorfield.dat');
  fid= fopen (name,'w');
  if (fid == -1)
    disp('error opening file '+name)
  end

  dlat = 1;
  for i = 1:nb_lat
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    nb_mesh = 3*nb_mesh;
    lat = 90-i*dlat;
    dlon = 360/(nb_mesh);
    for j=1:nb_mesh
      lon = j*dlon;
      count= fprintf (fid,'%f %f %f %f \n',lat,lon,u_s(i,j),v_s(i,j));
    end
  end
  fclose(fid);
end
