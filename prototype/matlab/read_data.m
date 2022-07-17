% Lit les concentrations sur la grille reduite
%

function [x] = read_data(k,nb_lat,nb_lat2)
  name = strcat('data/',num2str(k),'.dat')
  fid= fopen (name,'r');
  if (fid == -1)
    disp('error opening file ')
  end

  nb_cells = 3*get_nb_cells(nb_lat2,nb_lat,nb_lat2);
  precision = '16'; % Precision machine
  format = strcat('%.',precision,'*g %*g \n');

  tmp = fscanf (fid,'%g %*g %*g ');
  k = 1;
  x = zeros(nb_lat,nb_cells);
  for i=1:nb_lat
    % On garde seulement les valeurs au milieu de la maille
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    for j=1:nb_cells
      x(i,j) = tmp(k);
      k = k + 1;
    end
  end
  fclose(fid);
end
