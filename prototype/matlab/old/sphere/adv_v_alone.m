%
% advection nord-sud par secteur
%
function [q_s,f_v,mC_un_s]  =  adv_v_alone(q_s,duq_s,mCon_s,v_s,u_middle,u_step,...
    v_middle, v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2,dt);
  nb_sectors = 2;

  % Compute the meridional derivative
  function dq_s = compute_dq_s(q_s,duq_s,u_middle,u_step,nb_lat,nb_lat2,...
      nb_sectors,dlat)
    nb_mesh = 3*get_nb_mesh(nb_lat2,nb_lat,nb_lat2);
    dq_s = zeros(nb_lat,nb_mesh);
    for i = 2:nb_lat-1
      nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
      nb_mesh_prev = get_nb_mesh(i-1,nb_lat,nb_lat2);
      nb_mesh_next = get_nb_mesh(i+1,nb_lat,nb_lat2);

     % Up and down steps
      coslat1 = cos((90-dlat/2-(i-2)*dlat)*pi/180);
      dlon1 = 120*coslat1/(2*(i-1)-1);
      coslat2 = cos((90-dlat/2-i*dlat)*pi/180);
      dlon2 = 120*coslat2/(2*(i+1)-1);

      for j = 1:nb_mesh
        for nb_sec = 0:nb_sectors
          % Upper latitude
          coo1 = fix((u_middle(i,j))/u_step(i-1));
          dcoo1 = u_middle(i,j)-coo1*u_step(i-1);
          dcoo1 = dcoo1*dlon1;
          coo1 = coo1+1+nb_mesh_prev*nb_sec;
          q_sm = q_s(i-1,coo1)+(dcoo1-dlon1/2)*duq_s(i-1,coo1);

          % Lower latitude
          coo2 = fix((u_middle(i,j))/u_step(i+1));
          dcoo2 = u_middle(i,j)-coo2*u_step(i+1);
          dcoo2 = dcoo2*dlon2;
          coo2 = coo2+1+nb_mesh_next*nb_sec;
          q_sp = q_s(i+1,coo2)+(dcoo2-dlon2/2)*duq_s(i+1,coo2);

          %  Limitation de la pente 
          j2 = j+nb_mesh*nb_sec;
          dC_min = 2*min(abs(q_sp-q_s(i,j2)),abs(q_sm-q_s(i,j2)))/dlat;
          %
          dq_s(i,j2) = (q_sp-q_sm)*0.5/dlat;
          ypm = min(abs(dq_s(i,j2)),dC_min);
          dq_s(i,j2) = sign(dq_s(i,j2))*ypm;
      %    
        end
      end
    end
  end

  dlat = 1;
  %
  %  Limitation de la pente 
  dq_s = compute_dq_s(q_s,duq_s,u_middle,u_step,nb_lat,nb_lat2,nb_sectors,dlat);

  % Pour chaque noeud, on calcule l'advection entre jp (lat courante) et 
  % jn (lat suivante)
  nb_nodes = 3*get_nb_nodes(nb_lat2,nb_lat,nb_lat2);
  f_v = zeros(nb_lat,nb_nodes);
  for i = 1:nb_lat-1;
    coslat = cos((90-i*dlat)*pi/180);
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    nb_mesh_prev = get_nb_mesh(i,nb_lat,nb_lat2);
    nb_mesh_next = get_nb_mesh(i+1,nb_lat,nb_lat2);

    for j = 1:nb_nodes
      j_prev = i_mesh_prev(i,j);
      j_next = i_mesh_next(i,j);

      for nb_sec = 0:nb_sectors
        j2 = j + nb_sec*nb_nodes;
        jp = j_prev+nb_sec*nb_mesh_prev;
        ap = dq_s(i,jp);
        bp = q_s(i,jp)-ap*dlat*0.5;

        jm = j_next+nb_sec*nb_mesh_next;
        am = dq_s(i+1,jm);
        bm = q_s(i+1,jm)-am*dlat*0.5;
               
        v_val = v_s(i,j2);
        % Flux entre (i,jp) et (i+1,jm), associé à i+1
        % Si v < 0 : flux de (i+1,jm) vers (i,jp) (donc negatif)
        if (v_val < 0.)
          dy = -v_val*dt;
          f_v(i,j2) = (-am*dy*dy*0.5-bm*dy)*v_step(i,j)*120*coslat;
          f_v1(i,j2) = (-dy)*v_step(i,j)*120*coslat;;
        else
          dy = v_val*dt;
          f_v(i,j2) = (ap*0.5*(2*dlat-dy)*dy+bp*dy)*v_step(i,j)*120*coslat;; 
          f_v1(i,j2) = dy*v_step(i,j)*120*coslat;;
        end
      end
    end
  end

  nb_mesh = 3*get_nb_mesh(nb_lat2,nb_lat,nb_lat2);

  mCon_un_s = zeros(nb_lat,nb_mesh);
  for i=1:nb_lat
    for j=1:nb_mesh
      cur = q_s(i,j);
      if (cur ~= 0)
        mCon_un_s(i,j) = mCon_s(i,j)/cur;
      end
    end
  end
  
  for i = 1:nb_lat-1
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    nb_mesh_prev = get_nb_mesh(i,nb_lat,nb_lat2);
    nb_mesh_next = get_nb_mesh(i+1,nb_lat,nb_lat2);

    for j = 1:nb_nodes
      j_next = i_mesh_next(i,j);
      j_prev = i_mesh_prev(i,j);
      for nb_sec = 0:nb_sectors
        jm = j_next+nb_sec*nb_mesh_next;
        jp = j_prev+nb_sec*nb_mesh_prev;
       
        % Advection vers le sud
        j2 = j + nb_sec*nb_nodes;

        mCon_s(i+1,jm) = mCon_s(i+1,jm)+f_v(i,j2);
        mCon_s(i,jp) = mCon_s(i,jp)-f_v(i,j2);

        mCon_un_s(i+1,jm) = mCon_un_s(i+1,jm)+f_v1(i,j2);
        mCon_un_s(i,jp) = mCon_un_s(i,jp)-f_v1(i,j2);
      end
    end
  end

  for i = 1:nb_lat
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    n = 3*nb_mesh;
    for j = 1:n
      if (mCon_un_s(i,j) > 0)
        q_s(i,j) = mCon_s(i,j)/mCon_un_s(i,j);
      else
        q_s(i,j) = mCon_s(i,j);
      end
    end
  end
end
