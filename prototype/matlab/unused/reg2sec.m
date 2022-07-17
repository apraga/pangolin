% Interpolation d'une grille nb_lat x nb_lon vers une grille secteur

function [C_s,dC_s] = reg2sec(C,u_middle,u_step,nb_lat,nb_lat2)

  dlat = 90/nb_lat2;
  nb_lon = 3*(2*nb_lat2 - 1);

  for i=1:nb_lat
    % On garde seulement les valeurs au milieu de la maille
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    for j=1:nb_mesh
      pos = fix(u_middle(i,j)*120+0.5);
      pos = max(1,pos);
      for nb_sec=0:2
        C_s(i,nb_sec*nb_mesh + j) = C(i, nb_sec*120+ pos);
      end;
    end;
  end;
  % Calcul des derivées (2e ordre)
  % d C-s / d latitude
  dC_s(1,1)=0.;
  for i=2:nb_lat-1;
    nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
    for j=1:nb_mesh
     % Maille superieure et inférieure
      coo1 = fix(u_middle(i,j)/u_step(i-1)+0.5);
      coo1 = max(1,coo1);
      coo2 = fix(u_middle(i,j)/u_step(i+1)+0.5);
      coo2 = max(1,coo2);
      for nb_sec=0:2
        dC_s(i,j+nb_sec*120) = (C_s(i+1,coo2+nb_sec*120)-C_s(i-1,coo1+nb_sec*120))*0.5/dlat;    
      end;
    end;
  end;
  dC_s(nb_lat,1)=0.;
