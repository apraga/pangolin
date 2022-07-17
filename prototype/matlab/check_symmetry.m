% Check symmetry in respect to the equator line
function check_symmetry(k)
  nb_lat = 180;
  nb_lat2 = nb_lat/2;
  approx = 0;
  q_s = read_data(k,nb_lat,nb_lat2,approx);
  i_last = nb_lat;
  diff_q = -10000.;
  for i = 1:nb_lat2
    nb_mesh = 3*get_nb_mesh(i,nb_lat,nb_lat2,approx);
    for j = 1:nb_mesh
      diff_tmp = abs(q_s(i,j)-q_s(i_last,j));
      diff_q = max(diff_tmp,diff_q);
    end
    i_last = i_last - 1;
  end
  fprintf('Symmetry error %f \n',diff_q);
end


