% Correcting the v-component of the winds on the reduced grid
% We want mass preservation so the flux of v at a given latitude must be = 0

function [v_corr] = correction_v_red(v_s,v_step,nb_lat,nb_lat2)
  
  for i=1:nb_lat
    resv = 0.;
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);

    for j=1:nb_nodes
      step = v_step(i,j);
      for nb_sec=0:2;
        j2 = j + nb_sec*nb_nodes;
        resv = resv+v_s(i,j2)*step;
      end
    end
    % On divise par la surface totale = 3*1
    resv = resv/3;
    for j=1:3*nb_nodes
      v_corr(i,j) = v_s(i,j)-resv;
    end
  end
end
