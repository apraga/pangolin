%
% Correction de la composante u de la vitesse
% u_ss est le flux corrige, 
% div est la divergence de u (non corrige) sur une maille
% div2m est la divergence de u (non corrige) sur une latitude complete
function [u_ss,um,div,div2m] = cor_massu(u_s,v_ss,dcoo,nindm,nindp,nlat)
  %
  %
  dlat=90/nlat;
  pi=2*acos(0.);
  %
  for i=1:nlat;
    um(i)=0.;
    for j=1:2*i-1;
      div(i,j)=0.;
    end;
    for nsec=1:3;
      for j=1:2*i-1;
        j1=j+(nsec-1)*(2*i-1);
        um(i)=um(i)+u_s(i,j1);
        u_ss(i,j1)=0.;
      end;
    end;
  end;
  for n=1:3;
    div=div*0.;
    for i=1:nlat-1;
      coslat=cos((90-i*dlat)*pi/180);
      for j=1:4*i-1;
        vp=max(v_ss(i,j,n),0.);
        cd=vp;
        f_v1(i,j)=cd;
        vn=min(v_ss(i,j,n),0.);
        %
        cd=-vn;
        % flux de v :
        % la surface elementaire est rho cos(lat)d long
        % dons int(cos(lat) d long) = cos(lat) dcoo(i,j)
        f_v1(i,j)=(f_v1(i,j)-cd)*dcoo(i,j)*120*coslat;
      end;
    end;
    for j=1:4*nlat-1;
      f_v1(nlat,j)=0.;
    end;
    % div = somme des flux * surface
    for i=1:nlat-1;
      for j=1:4*i-1;
        jm=nindm(i,j);
        jp=nindp(i,j);
        div(i+1,jm)=div(i+1,jm)+f_v1(i,j);
        div(i,jp)=div(i,jp)-f_v1(i,j);
      end;
    end;
    for i=1:nlat;
      for j=1:2*i-1;
        jl=j+(n-1)*(2*i-1);
        div2(i,jl)=div(i,j);
      end;
    end;
    %
  end;
  % Calcul de u corrige
  %
  for i=1:nlat;
    div2m(i)=0.;
    u_ss(i,1)=u_s(i,1);
    for j=2:3*(2*i-1);
      u_ss(i,j)=u_ss(i,j-1)+div2(i,j-1)/dlat;
    end;
    for j=1:3*(2*i-1);
      div2m(i)=div2m(i)+div2(i,j);
    end;
  end;
