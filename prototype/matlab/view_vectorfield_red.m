% Displays winds at entire latitude (boundary of the cells, instead of the
% middle)
function view_vectorfield_red(v_middle,i_mesh_prev,nb_lat,nb_lat2,nb_lon,dlat,u_s,v_s)
  nb_max = 3*get_nb_nodes(nb_lat2-1,nb_lat,nb_lat2);
  n = nb_lat*nb_max;
  lat = zeros(n,1);
  lon = zeros(n,1);
  u = zeros(n,1);
  v = zeros(n,1);
  k = 1;
  for i=1:nb_lat-1
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
    for nb_sec=0:2
      for j=1:nb_nodes
        lat(k) = 90-i*dlat;
        lon(k) = (v_middle(i,j) + nb_sec)*120;
        %if (lon(k) > 180.)
        %  lon(k) = lon(k) - 360;
        %end
        j_prev = i_mesh_prev(i,j) + nb_cells*nb_sec;
        u(k) = u_s(i,j_prev);
        % V > 0 means south advection
        v(k) = -v_s(i,j + nb_nodes*nb_sec);
        k = k + 1;
      end
    end
  end
  quiver(lon,lat,u,v,2);
end
