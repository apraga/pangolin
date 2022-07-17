function n = get_nb_nodes(i,nb_lat,nb_lat2)
  n = 4*i - 1;
  if (i > nb_lat2)
    % Different du nb de mesh !
    i2 = nb_lat - i;
    n = 4*i2 - 1;
  end
end

