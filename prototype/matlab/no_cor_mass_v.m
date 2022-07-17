% Simple copie

function [v_s] = no_cor_mass_v(v,v_middle,v_step,nb_lat,nb_lat2,nb_lon)

  for i=1:nb_lat;
    resv = 0.;
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);

    for j=1:nb_nodes
      ind = fix((v_middle(i,j)-v_step(i,j)*0.5)*120)+1;
      for nb_sec=0:2;
        j2 = j + nb_sec*nb_nodes;
        ind2 = 120*nb_sec + ind;
        v_s(i,j2) = v(i,ind2);
      end
    end
end
