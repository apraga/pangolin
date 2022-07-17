% Calcul des coordonnees pour la grille reguliere
% middle contient le mileu des mailles
% step contient le pas entre les mailles
function [middle,step] = grille(nb_lat,nb_lat2)
  nb_lon = 2*nb_lat2 -1;
  step = zeros(nb_lat);
  middle = zeros(nb_lat,nb_lon);
  for i = 1:nb_lat
    nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
    % On se retreint a un seul secteur
    step(i) = 1./nb_cells;
    for j = 1:nb_cells
      middle(i,j) = (j-0.5)*step(i);
    end
  end
end
