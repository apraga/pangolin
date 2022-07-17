%
% Correction de la composante u de la vitesse
% u_corr est le flux corrige, 
function [u_corr] = correction_u(u_s,v_corr,v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2)

  nb_mesh = get_nb_mesh(nb_lat2,nb_lat,nb_lat2);
  div = zeros(nb_lat,3*nb_mesh);
  u_corr = zeros(nb_lat,3*nb_mesh);
  dlat = 180/nb_lat;
  nb_nodes = get_nb_nodes(nb_lat2,nb_lat,nb_lat2);
  f_v1 = zeros(nb_lat,3*nb_nodes);

  for i=1:nb_lat-1;
    coslat = cos((90-i*dlat)*pi/180);
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    for j=1:nb_nodes
      % flux de v :
      dS = v_step(i,j)*120*coslat;
      for nb_sec=0:2
        j2 = j + nb_nodes*nb_sec;
        f_v1(i,j2) = v_corr(i,j2)*dS;
      end
    end
  end

  % Somme des flux
  for i=1:nb_lat-1
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    nb_mesh_prev = get_nb_mesh(i,nb_lat,nb_lat2);
    nb_mesh_next = get_nb_mesh(i+1,nb_lat,nb_lat2);
    for j = 1:nb_nodes
      j_next = i_mesh_next(i,j);
      j_prev = i_mesh_prev(i,j);
      for nb_sec=0:2
        j_next2 = j_next + nb_sec*nb_mesh_next;
        j_prev2 = j_prev + nb_sec*nb_mesh_prev;
        j2 = j + nb_nodes*nb_sec;
        div(i+1,j_next2) = div(i+1,j_next2)+f_v1(i,j2);
        div(i,j_prev2) = div(i,j_prev2)-f_v1(i,j2);
      end 
    end
  end

  % Calcul de u corrige
  %
  %div(45,:)
  for i=1:nb_lat;
    u_corr(i,1) = u_s(i,1);
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    for j=2:3*nb_mesh
      u_corr(i,j) = u_corr(i,j-1)+div(i,j-1)/dlat;
    end
  end
