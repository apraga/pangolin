% Lit les concentrations sur la grille reduite
%

function [x] = ncl_read(k,nb_lat,nb_lat2)
  name = strcat('../../plot/data/',num2str(k),'.dat');
  fid= fopen (name,'r');
  if (fid == -1)
    disp('error opening file '+name)
  end

  nb_mesh = 3*get_nb_mesh(nb_lat2,nb_lat,nb_lat2);
  tmp = fscanf (fid,'%*g %*g %g ');
  k = 1;
  x = zeros(nb_lat,nb_mesh);
  %x = fscanf (fid,'%*g %*g %g ');
  for i=1:nb_lat
    % On garde seulement les valeurs au milieu de la maille
    nb_mesh = 3*get_nb_mesh(i,nb_lat,nb_lat2);
    for j=1:nb_mesh
      x(i,j) = tmp(k);
      k = k + 1;
    end
  end
  fclose(fid);
end
