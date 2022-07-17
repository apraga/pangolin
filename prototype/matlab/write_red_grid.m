function write_red_grid(v_middle,v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2)
name = strcat('grid_red.dat');
  fid= fopen (name,'w');
  %precision = '16'; % Precision machine
  %format = strcat('%.',precision,'g %.',precision,'g %.',precision,'g \n');
  format = strcat('%f %f %f %f \n');
  if (fid == -1)
    disp('error opening file ');
  end
  count= fprintf (fid,'prev next step middle \n');

  dlat = 1;

  for i = 1:nb_lat
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    for j=1:nb_nodes
      count= fprintf (fid,format,i_mesh_prev(i,j),i_mesh_next(i,j),v_step(i,j),v_middle(i,j));
    end
  end
  fclose(fid);

end
