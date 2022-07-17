% North-South advection
function [q_s]  =  adv_v_alone(q_s,duq_s,mCon_s,mCon_un_s,v_s,u_middle,u_step,...
    v_middle, v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2,dt,dlat);
  nb_sectors = 2;
  precision = eps*100;

  % Compute the meridional derivative
  function dq_s = compute_dq_s(q_s,duq_s,u_middle,u_step,nb_lat,nb_lat2,...
      nb_sectors,dlat)
    nb_cells = 3*get_nb_cells(nb_lat2,nb_lat,nb_lat2);
    dq_s = zeros(nb_lat,nb_cells);
    for i = 2:nb_lat-1
      nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
      nb_cells_prev = get_nb_cells(i-1,nb_lat,nb_lat2);
      nb_cells_next = get_nb_cells(i+1,nb_lat,nb_lat2);

     % Up and low steps
      coslat1 = cos((90-dlat/2-(i-2)*dlat)*pi/180);
      dlon1 = 120*coslat1/nb_cells_prev;
      coslat2 = cos((90-dlat/2-i*dlat)*pi/180);
      dlon2 = 120*coslat2/nb_cells_next;

      for j = 1:nb_cells
        u_middle_tmp(j) = u_middle(i,j);
      end
      for j = 1:nb_cells
        for nb_sec = 0:nb_sectors
          % Piecwise linear interpolation for the upper latitude
          u_mid = u_middle_tmp(j);
          coo1 = fix(u_mid/u_step(i-1));
          dcoo1 = u_mid-coo1*u_step(i-1);
          %dcoo1 = dcoo1*dlon1;
          dcoo1 = dcoo1*120*coslat1;
          coo1 = coo1+1+nb_cells_prev*nb_sec;
          q_up = q_s(i-1,coo1)+(dcoo1-0.5*dlon1)*duq_s(i-1,coo1);

          % Idem for lower
          coo2 = fix(u_mid/u_step(i+1));
          dcoo2 = u_mid-coo2*u_step(i+1);
          %dcoo2 = dcoo2*dlon2;
          dcoo2 = dcoo2*120*coslat2;
          coo2 = coo2+1+nb_cells_next*nb_sec;
          q_low = q_s(i+1,coo2)+(dcoo2-0.5*dlon2)*duq_s(i+1,coo2);

          %  Slope limitation (van Leer)
          j2 = j+nb_cells*nb_sec;
       %if (i == 81 && j2 == 483)
            %fprintf('for slope: %.16f %.16f\n', q_s(i-1,coo1), q_s(i+1, coo2));
%            fprintf('dq %.16f %.16f\n', duq_s(i-1,coo1), duq_s(i+1, coo2));
            %fprintf('dq: prev next %.16f %.16f\n', duq_s(i-1,coo1), duq_s(i+1,coo2));
            %fprintf('slope: prev next %.16f %.16f\n', q_up, q_low);
           % fprintf('diff prev next %.16f %.16f\n', dcoo1-0.5*dlon1, dcoo2-0.5*dlon2);
            %fprintf('qinterp %.16f %.16f\n', q_up, q_low);
          %end


          dq_s(i,j2) = slope_limit(q_up,q_s(i,j2),q_low,dlat);
          if (i == 79 && j2 == 117)
            %fprintf('slope limit %.16f', dq_s(i, j2));
          %%  fprintf('slope limit %15.10f prev cur next %15.10f %15.10f %15.10f\n', dq_s(i, j2), q_up, q_s(i,j2),q_low);
          end 
        end
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  dq_s = compute_dq_s(q_s,duq_s,u_middle,u_step,nb_lat,nb_lat2,nb_sectors,dlat);
  %fprintf('slope limit %15.10f\n', dq_s(100,1:2));

  % Pour chaque noeud, on calcule l'advection entre jp (lat courante) et 
  % jn (lat suivante)
  nb_nodes = 3*get_nb_nodes(nb_lat2,nb_lat,nb_lat2);
  % On s'arrete au pole nord
  for i = 1:nb_lat-1;
    coslat = cos((90-i*dlat)*pi/180);
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    nb_cells_prev = get_nb_cells(i,nb_lat,nb_lat2);
    nb_cells_next = get_nb_cells(i+1,nb_lat,nb_lat2);

    for j = 1:nb_nodes
      mesh_prev_tmp(j) = i_mesh_prev(i,j);
      mesh_next_tmp(j) = i_mesh_next(i,j);
      v_step_tmp(j) = v_step(i,j);
    end

    for nb_sec = 0:nb_sectors
      for j = 1:nb_nodes
        j_prev = mesh_prev_tmp(j);
        j_next = mesh_next_tmp(j);

        j2 = j + nb_sec*nb_nodes;
        jp = j_prev+nb_sec*nb_cells_prev;
        ap = dq_s(i,jp);
        bp = q_s(i,jp)-ap*dlat*0.5;

        jm = j_next+nb_sec*nb_cells_next;
        am = dq_s(i+1,jm);
        bm = q_s(i+1,jm)-am*dlat*0.5;

        v_val = v_s(i,j2);
        % Flux between (i,jp) and (i+1,jm), associated to i
        % if v < 0 : advection toward the north : from (i+1,jm) to (i,jp) ( so < 0)
        delta = v_step_tmp(j);
        if (v_val < 0.)
          dy = -v_val*dt;
          flux = (-am*dy*dy*0.5-bm*dy)*delta*120*coslat;
          flux_un = (-dy)*delta*120*coslat;
        else
          dy = v_val*dt;
          flux = (ap*0.5*(2*dlat-dy)*dy+bp*dy)*delta*120*coslat; 
          flux_un = dy*delta*120*coslat;
        end

        mCon_s(i+1,jm) = mCon_s(i+1,jm)+flux;
        mCon_s(i,jp) = mCon_s(i,jp)-flux;

        mCon_un_s(i+1,jm) = mCon_un_s(i+1,jm)+flux_un;
        mCon_un_s(i,jp) = mCon_un_s(i,jp)-flux_un;

      end
    end
  end
  %fprintf('mass %.16f %.16f \n', mCon_s(79,117),mCon_un_s(79,117));
  for i = 1:nb_lat
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    for j = 1:nb_cells
      test = abs(mCon_un_s(i,j)) > precision;
      if (test)
        q_s(i,j) = mCon_s(i,j)/mCon_un_s(i,j);
      else
        q_s(i,j) = 0;
      end
    end
  end
end
