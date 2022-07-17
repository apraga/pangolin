%
% calcul des coordonnees des noeuds
% co3m contient le milieu de 2 noeuds
% dco3 contient la coordonnees d'un noeuds en fonction du precedent
% indp contient le numero de la maille superieure a laquelle le noeud appartient
% indm contient le numero de la maille inferieure a laquelle le noeud appartient
function [co3m,dco3,nsec,indp,indm] = sub_grille(nlat)
  %
  %nlat=90;
  %
  nsec=4*nlat-1;
  for n=1:nlat;
    for ia=1:2*n-1;
      co1(n,ia)=(ia-1)/(2*n-1);
    end;
    for ib=1:2*n+1;
      co2(n,ib)=(ib-1)/(2*n+1);
    end;
  end;
  xco1=co1;
  xco2=co2;
  for n=1:nlat;
    co3(n,1)=0.;
    for ic=2:4*n-1;
      mco1=1.;
      mco2=1.;
      for ia=2:2*n-1;
        mco1=min(mco1,co1(n,ia));
      end;
      for ib=2:2*n+1;
        mco2=min(mco2,co2(n,ib));
      end;
      mm=min(mco1,mco2);
      co3(n,ic)=co3(n,ic-1)+mm;
      for ia=1:2*n-1;
        if(co1(n,ia) < 0.99999999999);
          co1(n,ia)=co1(n,ia)-mm;
        end;
        if(co1(n,ia) < 1e-20);
          co1(n,ia)=1.;
        end;
      end;
      for ib=1:2*n+1;
        if(co2(n,ib) < 0.99999999999);
          co2(n,ib)=co2(n,ib)-mm;
        end;
        if(co2(n,ib) < 1e-20);
          co2(n,ib)=1.;
        end;
      end;
    end;
    co3(n,4*n)=1.;
  end;
  minco3=1.;
  for n=1:nlat;
    for ic=2:4*n;
      dco3(n,ic-1)=co3(n,ic)-co3(n,ic-1);
      co3m(n,ic-1)=(co3(n,ic)+co3(n,ic-1))*0.5;
      minco3=min(minco3,dco3(n,ic-1));
    end;
  end;
  %
  for n=1:nlat;
    for ic=1:4*n-1;
      indp(n,ic)=fix(co3m(n,ic)*(2*n-1))+1; % Va de 1 a 2n-1
      indm(n,ic)=fix(co3m(n,ic)*(2*n+1))+1;
    end;
  end;
