% Advection sur toute la sphère

nb_lat = 180;
nb_lat2 = nb_lat/2;
nb_lon = 360;
% Rayon de la planete
radius = 1.;
% Conditions initiales
% distrib = 1 : advection zonale
% distrib = 2 : advection meridionale
% distrib = 3 : advection bi_dim (Hourdin)
% distrib = 4 : advection bi_dim
distrib = 3; 

if (distrib == 1)
  dt = 120;
elseif (distrib == 3)
  dt = 100;
else
  dt = 20;
end
q = init_C(nb_lat,nb_lat2,nb_lon,distrib);

[u_middle,u_step] = grille(nb_lat,nb_lat2);

[v_middle,v_step,nb_nodes,i_mesh_prev,i_mesh_next] = sub_grille(nb_lat,nb_lat2);

% Initialisation de la circulation
[u,v] = circul(nb_lat,nb_lon,radius,distrib);
[u_s] = reg2sec_u(u,u_middle,u_step,nb_lat,nb_lat2);

%[v_corr] = no_cor_mass_v(v,v_middle,v_step,nb_lat,nb_lat2,nb_lon);
% Correction de u et v pour conservation de la masse
[v_corr] = correction_v(v,v_middle,v_step,nb_lat,nb_lat2,nb_lon);
[u_corr] = correction_u(u_s,v_corr,v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2);
v_s = v_corr;
u_s = u_corr;
print_maxmin(u_s,v_s,'corr');
%write_vectorfield(u_s,v_s,nb_lat,nb_lat2,nb_lon,v_step,v_middle,i_mesh_prev,i_mesh_next);


[q_s,dq_s]= reg2sec(q,u_middle,u_step,nb_lat,nb_lat2);
ncl_write(0,q_s,nb_lat,nb_lat2,nb_lon);

% Calcul de la masse
[m0] = masse(q_s,nb_lat,nb_lat2,nb_lon);

%% Integration temporelle 
%[q_s] = ncl_read(300,nb_lat,nb_lat2);
error = 0.1;
for k=1:1500
  % Advection est-ouest de u
  [q_s,duq_s,mCon_s] = adv_u_alone(q_s,u_s,nb_lat,nb_lat2,dt);

  % Advection nord-sud de v
  [q_s,f_v] = adv_v_alone(q_s,duq_s,mCon_s,v_s,u_middle,u_step,v_middle,...
    v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2,dt);

  k
  if (mod(k,50) == 0)
    ncl_write(k,q_s,nb_lat,nb_lat2,nb_lon);
    [m] = masse(q_s,nb_lat,nb_lat2);
    if (abs(m - m0) > error)
      disp('divergence !');
      return
    end
  end
end

% diagnostic masse
[m0] = masse(q_s,nb_lat,nb_lat2);
