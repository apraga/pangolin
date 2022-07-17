function write_vectorfield(u_s,v_s,nb_lat,nb_lat2,nb_lon,v_step,i_mesh_prev,...
    i_mesh_next)
  name = strcat('../../plot/data/data_vectorfield.dat');
  fid= fopen (name,'w');
  if (fid == -1)
    disp('error opening file '+name)
  end

  dlat = 1;
  for i = 2:nb_lat
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    lat = 90-i*dlat;
%    dlon = 360/(nb_mesh);
    sum_lon = 0;
    for j=1:nb_nodes
      j_prev = i_mesh_prev(i,j);
      j_next = i_mesh_next(i,j);
      for nb_sec=0:2
        sum_lon = sum_lon + v_step(i,j);
        lon = sum_lon*120 + nb_sec*120;
        j_u = j_prev + nb_sec*nb_mesh;
        j_v = j + nb_sec*nb_nodes;
        count= fprintf (fid,'%f %f %f %f \n',lat,lon,u_s(i,j_u)*5,v_s(i,j_v)*5);
      end
    end
  end
  fclose(fid);
end
