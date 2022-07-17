% Nubmer of cells at a latitude i for a sector
% Linear approximation ( n = 2i-1)
function n = get_nb_cells(i,nb_lat,nb_lat2)
  i2 = i;
  if (i > nb_lat2)
    i2 = nb_lat + 1 - i;
  end

  n = 2*i2 - 1;
end


