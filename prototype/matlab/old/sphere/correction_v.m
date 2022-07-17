% Correction de la composante v de la vitesse:
% flux a latitude constante = 0

function [v_corr] = correction_v(v,v_middle,v_step,nb_lat,nb_lat2,nb_lon)
  %
  for i=1:nb_lat;
    resv = 0.;
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);

    for j=1:nb_nodes
      ind = fix((v_middle(i,j)-v_step(i,j)*0.5)*120);
      step = v_step(i,j);
      for nb_sec=0:2;
        j2 = j + nb_sec*nb_nodes;
        ind2 = 120*nb_sec + max(1,ind);
        v_corr(i,j2) = v(i+1,ind2);
        resv = resv+v_corr(i,j2)*step;
      end
    end
    % On divise par la surface totale = 3*1
    resv = resv/3;
    for j=1:nb_nodes
      for nb_sec=0:2;
        j2 = j + nb_sec*nb_nodes;
        v_corr(i,j2) = v_corr(i,j2)-resv;
      end
    end
  end
