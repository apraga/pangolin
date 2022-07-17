%
% calcul des coordonnees des noeuds
% middle contient le milieu de 2 noeuds
% step contient la coordonnees d'un noeuds en fonction du precedent
% mesh_prev contient le numero de la maille superieure a laquelle le noeud appartient
% mesh_next contient le numero de la maille inferieure a laquelle le noeud appartient
function [middle,step,nb_nodes,mesh_prev,mesh_next] = sub_grille(nb_lat,nb_lat2)
  nb_nodes = get_nb_nodes(nb_lat2,nb_lat,nb_lat2);

  mesh_prev = zeros(nb_lat,nb_nodes);
  mesh_next = zeros(nb_lat,nb_nodes);
  mesh_next = zeros(nb_lat,nb_nodes);
  step = zeros(nb_lat,nb_nodes);
  middle = zeros(nb_lat,nb_nodes);

  for n=1:nb_lat2;
    for ia=1:2*n-1;
      co1(n,ia) = (ia-1)/(2*n-1);
    end;
    for ib=1:2*n+1;
      co2(n,ib) = (ib-1)/(2*n+1);
    end;
  end;
  for n=1:nb_lat2;
    co3(n,1) = 0.;
    for ic = 2:4*n-1;
      mco1 = 1.;
      mco2 = 1.;
      for ia=2:2*n-1;
        mco1 = min(mco1,co1(n,ia));
      end;
      for ib=2:2*n+1;
        mco2 = min(mco2,co2(n,ib));
      end;
      mm = min(mco1,mco2);
      co3(n,ic) = co3(n,ic-1)+mm;
      for ia=1:2*n-1;
        if(co1(n,ia) < 0.99999999999);
          co1(n,ia) = co1(n,ia)-mm;
        end;
        if(co1(n,ia) < 1e-20);
          co1(n,ia) = 1.;
        end;
      end;
      for ib=1:2*n+1;
        if(co2(n,ib) < 0.99999999999);
          co2(n,ib) = co2(n,ib)-mm;
        end;
        if(co2(n,ib) < 1e-20);
          co2(n,ib) = 1.;
        end;
      end;
    end;
    co3(n,4*n) = 1.;
  end;
  for n=1:nb_lat2;
    for ic = 2:4*n;
      step(n,ic-1) = co3(n,ic)-co3(n,ic-1);
      middle(n,ic-1) = (co3(n,ic)+co3(n,ic-1))*0.5;
    end;
  end;
  
  for n=1:nb_lat2;
    for ic=1:4*n-1;
      mesh_prev(n,ic) = fix(middle(n,ic)*(2*n-1))+1; % Va de 1 a 2n-1
      mesh_next(n,ic) = fix(middle(n,ic)*(2*n+1))+1;
    end;
  end;

  % Sur l'equateur, on connecte les 2 hemisphères
  mesh_next(nb_lat2,:) = mesh_prev(nb_lat2,:);

  % Hemisphere sud
  i_prev = nb_lat2-1;
  for i = nb_lat2+1:nb_lat-1
    mesh_prev(i,:) = mesh_next(i_prev,:);
    mesh_next(i,:) = mesh_prev(i_prev,:);
    step(i,:) = step(i_prev,:);
    middle(i,:) = middle(i_prev,:);
    i_prev = i_prev - 1;
  end
  %mesh_prev(nb_lat2,:)
  %mesh_prev(nb_lat,:) = mesh_next(i_prev,:);
  %mesh_next(nb_lat,:) = mesh_prev(i_prev,:);

end
