%C_liste vers C_secteur
function [C_s] = Cl2Cs(C_l,nb_lat,nsec);
  %
  sum_mesh = 0;
  for i = 1:nb_lat
    nb_mesh = get_nb_mesh(i,nb_lat);
    n = 3*nb_mesh;
    for j = 1:nb_mesh
      for nb_sec = 0:2
        n = nb_sec*nb_mesh+j;
        C_s(i,n) = C_l(sum_mesh+n);
      end
    end
    sum_mesh = sum_mesh + n;
  end
end
