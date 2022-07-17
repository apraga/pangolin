%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advection on the sphere with a Van Leer scheme on an area-preserving grid
% Main file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We deduce the number of latitude
nb_lat2 = 80;
nb_lat = 2*nb_lat2;

% Init distribution and winds circulation according to the situation
% {meridional,zonal,bidim_hourdin,solid_rotation,
% standard_rough}
scenario = 'standard_cosine';
%
% Resume the computation
k_restart = 0;
k_end = 3;%1200; % Default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correct winds for mass conservation (default)
correction = 1;


fprintf('Scenario : %s\n',scenario);
switch scenario
  case 'zonal'
    distrib = 2; 
    winds = 1; 
    correction = 0;
  case 'meridional'
    distrib = 1; 
    winds = 2; 
    correction = 0;
  case 'bidim_hourdin'
    winds = 3; 
    distrib = 8; 
  case 'solid_rotation'
    distrib = 7; 
    winds = 4; 
  case 'standard_rough'
    distrib = 5; 
    winds = 5; 
  case 'standard_gaussian'
    distrib = 3; 
    winds = 5; 
  case 'standard_dv'
    distrib = 6; 
    winds = 6; 
  case 'standard_cosine'
    distrib = 4; 
    winds = 5; 
  case 'standard_cosine_corr'
    distrib = 6; 
    winds = 5; 
  otherwise
    distrib = 1; 
    winds = 1; 
end

% Constant delta latitude
dlat = 90./nb_lat2;

[u_middle,u_step] = grille(nb_lat,nb_lat2);

[v_middle,v_step,nb_nodes,i_mesh_prev,i_mesh_next] = sub_grille(nb_lat,...
  nb_lat2);

%write_red_grid(v_middle,v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2)
% Wind initialisation
t = 0.;
period = 5.; % For standart tests
[u_s,v_s,dt] = circul_red(t,v_middle,nb_lat,nb_lat2,winds,dlat,period);

%view_vectorfield_red(v_middle,i_mesh_prev,nb_lat,nb_lat2,dlat,u_s,v_s);

%write_winds(u_s, v_s, v_middle, nb_lat, nb_lat2, dlat);
if (correction)
  [v_s] = correction_v_red(v_s,v_step,nb_lat,nb_lat2);
  [u_s] = correction_u_red(u_s,v_s,v_step,i_mesh_prev,i_mesh_next,nb_lat,...
    nb_lat2,dlat);
end
%write_winds(u_s, v_s, v_middle, nb_lat, nb_lat2, dlat);


% Init on reduced grid
q_s = init_C_red(nb_lat,nb_lat2,distrib,dlat);
% Total mass
m0 = masse(q_s,nb_lat,nb_lat2,0,0,1,dlat);

[C,k_end] = estim_C(u_s,v_s,dt,nb_lat,nb_lat2,k_end);

if (strfind(scenario,'standard'))
  k_end = period/dt;
  fprintf('Overriding kend : new value =%f\n', k_end);
end
fprintf(' dt=%.16f\n', dt);

%k_end = 1250;

if (k_restart < 1)
  write_data(0,q_s,nb_lat,nb_lat2,dlat);
  if (C > 1+ eps)
    fprintf(' Courant number > 1\n');
  end
else
  % Resuming previous computation
  q_s = read_data(k_restart,nb_lat,nb_lat2);
  m0 = masse(q_s,nb_lat,nb_lat2,m0,0,0,dlat);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Iteration from %d to %d, dt=%d \n',k_restart,k_end, dt);

for k=k_restart+1:k_end
  fprintf('-----------------iteration %d----------------- \n',k);
  % East - west advection
  %fprintf('before zonal %.16f\n', q_s(80,118));
  %fprintf('before zonal %.16f\n', q_s(79,117));
  [q_s,duq_s,mCon_s,mCon_un_s] = adv_u_alone(q_s,u_s,nb_lat,nb_lat2,dt,dlat);
  
  %fprintf('before merid %.16f %.16f\n', q_s(80,1), q_s(80,477));

  % North-South advection
  q_s = adv_v_alone(q_s,duq_s,mCon_s,mCon_un_s,v_s,u_middle,u_step,v_middle,...
    v_step,i_mesh_prev,i_mesh_next,nb_lat,nb_lat2,dt,dlat);

  if (mod(k,50) == 0)
    write_data(k,q_s,nb_lat,nb_lat2,dlat);
    m = masse(q_s,nb_lat,nb_lat2,m0,0,1,dlat);
  end
  t = double(k*dt);
  [u_s,v_s] = update_winds(u_s,v_s,t,dt,v_step,v_middle,i_mesh_prev,i_mesh_next,...
    nb_lat,nb_lat2,winds,correction,dlat);
  %fprintf('%d %d \n',q_s(41,1),q_s(41,2))
end

% Position of the max
%[i,j] = ind2sub(size(q_s), find(q_s>1))
fprintf('max %.16f %d %d\n', max(max(q_s)))
fprintf('min %.16f\n',min(min(q_s)))

write_data(k,q_s,nb_lat,nb_lat2,dlat);
m = masse(q_s,nb_lat,nb_lat2,m0,0,1,dlat);
