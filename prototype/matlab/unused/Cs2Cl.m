%C_secteur vers C_liste
% On stocke d'abord a latitude constante
function [C_l,xlat,xlon] = Cs2Cl(C_s,nb_lat);
  % xlat et xlon sont les milieux des mailles 
  % n^2 mailles * 3 secteurs
  nb_lat2 = (nb_lat+1)/2;
  n2 = 6*nb_lat2*nb_lat2;

  C_l = zeros(1,n2);
  xlat = zeros(1,n2);
  xlon = zeros(1,n2);
  %
  dlat = 179/nb_lat;
  sum_mesh = 0;
  for i=1:nb_lat;
    nb_mesh = get_nb_mesh(i,nb_lat);
    dlon = 360/nb_mesh;
    %  n = 3*(i-1)*(i-1);
    for j=1:nb_mesh
      % tailles de toutes les mailles stockees jusqu'a la latitude i
      for nb_sec = 0:2
        n = nb_sec*nb_mesh+j;
        C_l(sum_mesh+n) = C_s(i,n);
        xlat(sum_mesh+n) = 90-dlat/2-dlat*(i-1);
        xlon(sum_mesh+n) = dlon/2+dlon*(n-1);
      end
    end
    sum_mesh = sum_mesh + 3*nb_mesh;
  end
