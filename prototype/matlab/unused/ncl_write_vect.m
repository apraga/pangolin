function ncl_write_vect(name,u,v,nb_lat,nb_lon)
  fid= fopen (name,'w');
  if (fid == -1)
    disp('error opening file '+name)
  end

  nb_lat2 = (nb_lat+1)/2;
  dlat = 1;

  for i = 1:nb_lat
    nb_mesh = 3*(2*i - 1);
    if (i > nb_lat2)
      nb_mesh = 3*(2*(nb_lat+1-i)-1);
    end
    lat = 90-i*dlat;
    dlon = 360/(nb_mesh);
    for j=1:nb_mesh
      lon = j*dlon;
      count= fprintf (fid,'%f %f %f %f \n',lat,lon,u(i,j),v(i,j));
    end
  end
  fclose(fid);
end
