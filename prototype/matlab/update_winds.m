% Update winds at the next timestep
function [u_s,v_s] = update_winds(u_s,v_s,t,dt,v_step,v_middle,i_mesh_prev,i_mesh_next,...
    nb_lat,nb_lat2,winds,correction,dlat)
  can_correct = false;
  if(winds == 5)
    period = 5.;
    [u_s,v_s] = winds_standard_cv(u_s,v_s,t,dt,period,nb_lat,nb_lat2,v_middle,dlat);
    can_correct = true;
  elseif (winds == 6)
    period = 5.;
    [u_s,v_s] = winds_standard_dv(u_s,v_s,t,dt,period,nb_lat,nb_lat2,v_middle,dlat);
    can_correct = true;
  end

  if(can_correct && correction)
    [v_s] = correction_v_red(v_s,v_step,nb_lat,nb_lat2);
    [u_s] = correction_u_red(u_s,v_s,v_step,i_mesh_prev,i_mesh_next,nb_lat,...
      nb_lat2,dlat);
  end
end
