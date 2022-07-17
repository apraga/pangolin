% Simple copie

function [v_s] = no_cor_mass_v(v,v_middle,v_step,nb_lat,nb_lat2,nb_lon)

  %v_s = zeros(nb_lat,nb_lon)
  for i=1:nb_lat
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    for j=1:nb_nodes
      extr1 = fix((v_middle(i,j)-v_step(i,j)*0.5)*120);
      for nb_sec=0:2
        extr = 120*nb_sec+max(1,extr1);
        j1 = j+nb_sec*nb_nodes;
        v_s(i,j1) = v(i,extr);
      end;
    end;
  end
end
