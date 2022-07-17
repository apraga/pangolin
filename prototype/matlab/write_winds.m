% Write winds without interpolation  on the reduced grid
function write_winds(u, v, v_middle, nb_lat, nb_lat2, dlat)
  function fid = open_data(name)
    fid= fopen (name,'w');
    if (fid == -1)
      disp('error opening file ');
    end
  end

  format = '%20.13f %20.13f %20.13f\n';
  fid = open_data('data/u_corr.dat');
  % Write u
  for i = 1:nb_lat
    nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
    nb_cells = 3*nb_cells;
    lat = 90-(i-0.5)*dlat;
    dlon = 360/(nb_cells);
    for j=1:nb_cells
      lon = (j-1)*dlon;
      count= fprintf (fid, format, u(i, j), lat, lon);
    end
    % Periodic conditions
    lon = nb_cells*dlon;
    count= fprintf (fid, format, u(i, 1), lat, lon);
 
  end
  fclose(fid);

  fid = open_data('data/v_corr.dat');
  % Write v
  for i = 1:nb_lat
    nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
    nb_cells = 3*nb_cells;
    lat = 90.-i*dlat;
    dlon = 360/(nb_cells);
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    for nb_sec = 0:2
      for j=1:nb_nodes
        lon = v_middle(i, j)*120. + 120.*nb_sec;
        j2 = j + nb_nodes*nb_sec;
        count= fprintf (fid, format, v(i, j2), lat, lon);
      end
    end 
  end

  % Write v
  fclose(fid);
end
