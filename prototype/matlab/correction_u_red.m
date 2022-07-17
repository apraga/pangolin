% v being already corrected, we correct u on the reduced grid
% We want the divergence of the winds = 0

function [u_corr] = correction_u_red(u_s,v_corr,v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2,dlat)

  nb_cells = get_nb_cells(nb_lat2,nb_lat,nb_lat2);
  div = zeros(nb_lat,3*nb_cells);
  u_corr = zeros(nb_lat,3*nb_cells);
  nb_nodes = get_nb_nodes(nb_lat2,nb_lat,nb_lat2);
  f_v1 = zeros(nb_lat,3*nb_nodes);

  % Somme des flux
  for i=1:nb_lat-1
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    nb_cells_prev = get_nb_cells(i,nb_lat,nb_lat2);
    nb_cells_next = get_nb_cells(i+1,nb_lat,nb_lat2);
    coslat = cos((90-i*dlat)*pi/180);
    for j = 1:nb_nodes
      j_next = i_mesh_next(i,j);
      j_prev = i_mesh_prev(i,j);
      dS = v_step(i,j)*120*coslat;
      for nb_sec=0:2
        j2 = j + nb_nodes*nb_sec;
        flux = v_corr(i,j2)*dS;

        j_next2 = j_next + nb_sec*nb_cells_next;
        j_prev2 = j_prev + nb_sec*nb_cells_prev;
        div(i+1,j_next2) = div(i+1,j_next2) + flux;
        div(i,j_prev2) = div(i,j_prev2) - flux;
      end 
    end
  end

  % Calcul de u corrige
  for i=1:nb_lat;
    u_corr(i,1) = u_s(i,1);
    nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
    for j=2:3*nb_cells
      u_corr(i,j) = u_corr(i,j-1)+div(i,j-1)/dlat;
    end
  end
end
