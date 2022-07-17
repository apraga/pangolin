%
% Interpolation d'une grille nb_lat x nb_lon+1 vers une grille secteur
%
function [u_s] = reg2sec_u(u,u_middle,u_step,nb_lat,nb_lat2)
  % u est defini seulement a l'interface
  nb_lon = 2*nb_lat2-1;
  u_s = zeros(nb_lat,nb_lon);

  for i=1:nb_lat
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    for j=1:nb_mesh
      % extr est une extremite de la maille (en degre) sur une latitude donnee
      extr = fix((u_middle(i,j)-u_step(i)*0.5)*120);
      for nb_sec=0:2;
        extr2 = 120*nb_sec+max(1,extr);
        j1 = j+nb_sec*nb_mesh;
        u_s(i,j1) = u(i,extr2);
      end
    end
  end
end
