function n = get_nb_nodes(i,nb_lat,nb_lat2)
  i2 = i;

  % Connection between hemispheres at the equator
  if (i == nb_lat2)
    n = get_nb_cells(i,nb_lat,nb_lat2);
    return;
  elseif (i > nb_lat2)
    i2 = nb_lat - i;
  end
  
  n = 4*i2 - 1;
end

