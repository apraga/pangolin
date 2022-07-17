% Calcul de la masse
%
function [mass] = masse(C_s,nb_lat,nb_lat2,nb_lon)
  %
  % Calcul des flux

  dlat = 90/nb_lat2;
  mass = 0;
  for i=1:nb_lat;
    coslat = cos((90-dlat/2-(i-1)*dlat)*pi/180);
    nb_mesh = 3*get_nb_mesh(i,nb_lat,nb_lat2);
    dlon = 360*coslat/nb_mesh;
    dS = dlon*dlat;
    for j=1:nb_mesh;
      % concentration = masse / surface 
      mass = mass+C_s(i,j)*dlat*dlon;
    end
  end
  fprintf('mass = %.4f \n',mass);
end
