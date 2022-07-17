% Courant number estimation 
function [max_C,k_end_new] = estim_C(u_s,v_s,dt,nb_lat,nb_lat2,nb_lon,k_end)
  %u_max = max(max(u_s));
  %u_min = min(min(u_s));
  %v_max = max(max(v_s));
  %v_min = min(min(v_s));

  dy = 1.;
  nb_cells = 3*get_nb_cells(nb_lat2,nb_lat,nb_lat2);
  estim = zeros(nb_lat2,nb_cells);
  max_u = -10000;
  max_v = -10000;

  for i=1:nb_lat2
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    dx = 360./nb_cells;
    for j=1:nb_cells
      max_u = max(max_u,abs(u_s(i,j))/dx);
    end
    nb_nodes = 3*get_nb_cells(i,nb_lat,nb_lat2);
    for j=1:nb_nodes
      max_v = max(max_v,abs(v_s(i,j))/dy);
    end
  end
  max_C = max(max_u*dt,max_v*dt);

  fprintf('estimation C %f\n',max_C);
  maxu = max(max(u_s));
  if (maxu > 0)
    k_end_new = int32(360./(maxu*dt));
%    fprintf('estimation period %d\n',k_end_new);
    fprintf(' New k_end : %d\n',k_end_new);
  else
    k_end_new = k_end;
  end

end
