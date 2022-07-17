% Check winds symmetry
function check_winds(u_s,v_s,nb_lat,nb_lat2)
  i2 = nb_lat-1;
  diff_v = -10;
  for i = 1:nb_lat2-1
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    for j = 1:nb_cells
      diff_tmp_v = abs(v_s(i,j)-v_s(i2,j));
      diff_v = max(diff_tmp_v,diff_v);
    end
    i2 = i2 - 1;
  end

  i2 = nb_lat;
  diff_u = -10;
  for i = 1:nb_lat2
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    for j = 1:nb_cells
      diff_tmp_u = abs(u_s(i,j)+u_s(i2,j));
      diff_u = max(diff_tmp_u,diff_u);
    end
    i2 = i2 - 1;
  end
  fprintf('Symmetry error v %f , u %f \n',diff_v,diff_u);
end


