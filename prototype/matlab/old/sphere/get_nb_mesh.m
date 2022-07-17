% Nombre de mesh pour une latitude i pour un secteur

function n = get_nb_mesh(i,nb_lat,nb_lat2)
  n = 2*i - 1;
  if (i > nb_lat2)
    i2 = nb_lat + 1 - i;
    n = 2*i2 - 1;
  end
end


