%
% advection nord-sud par secteur
% mC_l : masse du traceur 
% mC_u : masse unitaire du traceur
% duC_l : derivee de la concentration
function [mC_l,mC_un,duC_l] = adv_u_alone(C_l,u_s,mco,dco,nlat,nlon,dt)
  %
  %
  % Calcul des flux
  %
  dlat=90/nlat;
  pi=2*acos(0.);
  mass=0;
  for i=1:nlat;
    %for i=2:nlat;
    coslat=cos((90-dlat/2-(i-1)*dlat)*pi/180);
    jl=3*(2*i-1);
    for j=1:jl;
      ind=3*(i-1)*(i-1)+j;
      Con(j)=C_l(ind);
      xun(j)=1;
    end;
    %
    % Calcul des dérivées
    %
    dlon=360*coslat/jl;
    for j=2:jl-1;
      % Limitation de la pente (van Leer)
      dCon(j)=(Con(j+1)-Con(j-1))*0.5/dlon;
      dC_min=2*min(abs(Con(j+1)-Con(j)),abs(Con(j-1)-Con(j)))/dlon;
      ypm=min(abs(dCon(j)),dC_min);
      extr=1.;
      dCon(j)=sign(dCon(j))*ypm*extr;
    end;
    % Conditions aux limites périodiques
    dCon(1)=(Con(2)-Con(jl))*0.5/dlon;
    dC_min=2*min(abs(Con(2)-Con(1)),abs(Con(jl)-Con(1)))/dlon;
    ypm=min(abs(dCon(1)),dC_min);
    dCon(1)=sign(dCon(1))*ypm*extr;
    %
    dCon(jl)=(Con(1)-Con(jl-1))*0.5/dlon;
    dC_min=2*min(abs(Con(1)-Con(jl)),abs(Con(jl-1)-Con(jl)))/dlon;
    ypm=min(abs(dCon(jl)),dC_min);
    dCon(jl)=sign(dCon(jl))*ypm*extr;
    %
    % Attention les concentrations sont stockees aux milieu des mailles
    for j=2:jl;
      am=dCon(j-1);
      % Concentration en j-1 en fait
      bm=Con(j-1)-am*0.5*dlon;
      ap=dCon(j);
      % Concentration en j
      bp=Con(j)-ap*0.5*dlon;
      up=max(u_s(i,j),0.);
      cd=up*dt;
      % le signe de u impose de sens d'advection
      % si u > 0, q chapeau = q_{j-1} + am*(dlon - 1/2 cd)
      f_u(i,j)=am*0.5*(2*dlon-cd)*cd+bm*cd; 
      f_u1(i,j)=cd;
      un=min(u_s(i,j),0.);
      cd=-un*dt;
      % si u < 0, q chapeau = q_j - ap*(1/2 cd)
      f_u(i,j)=f_u(i,j)-ap*cd*cd*0.5-bp*cd;
      f_u1(i,j)=f_u1(i,j)-cd;
    end;
    % Traitement d'un bord
    j=1;
    am=dCon(jl);
    bm=Con(jl)-am*0.5*dlon;
    ap=dCon(j);
    bp=Con(j)-ap*0.5*dlon;
    up=max(u_s(i,j),0.);
    cd=up*dt;
    f_u(i,j)=am*0.5*(2*dlon-cd)*cd+bm*cd; % Pourquoi 2 ?
    f_u1(i,j)=cd;
    un=min(u_s(i,j),0.);
    cd=-un*dt;
    f_u(i,j)=f_u(i,j)-ap*cd*cd*0.5-bp*cd;
    f_u1(i,j)=f_u1(i,j)-cd;
    %
   
    for j=1:jl-1;
      % Variation de masse
      dx(j)=(f_u(i,j+1)-f_u(i,j))*dlat;
      dxun(j)=(f_u1(i,j+1)-f_u1(i,j))*dlat;
      % Mise à jour des masses : on enlève la masse advecté (si u > 0)
      % sinon on l'ajoute
      mCon(j)=-dx(j)+Con(j)*dlat*dlon;
      mxun(j)=-dxun(j)+xun(j)*dlat*dlon;
      %
    end;
    % Traitement de l'autre bord
    j=jl;
    dx(j)=(f_u(i,1)-f_u(i,j))*dlat;
    dxun(j)=(f_u1(i,1)-f_u1(i,j))*dlat;
    mCon(j)=-dx(j)+Con(j)*dlat*dlon;
    mxun(j)=-dxun(j)+xun(j)*dlat*dlon;
    %if (i ==1)
    %  %dCon(1:3)
    %  fprintf('Flux %f %f %f \n',f_u(i,1),f_u(i,2),f_u(i,3));
    %  fprintf('mCo  %f %f %f \n',mCon(1),mCon(2),mCon(3));
    %  fprintf('dx  %f %f %f \n',dx(1),dx(2),dx(3));
    %   fprintf('dlat dlon %f %f \n',dlat,dlon);

    %end
 
    %
    % Calcul des dérivées
    % On reorganise la masse dans les cellules
    dS = dlon/dlat;
    for j=1:jl;
      Con(j)=mCon(j)/dS;
    end;
    %
    for j=2:jl-1;
      dCon(j)=(Con(j+1)-Con(j-1))*0.5/dlon;
    end;
    dCon(1)=(Con(2)-Con(jl))*0.5/dlon;
    dCon(jl)=(Con(1)-Con(jl-1))*0.5/dlon;
    %
    for j=1:jl;
      ind=3*(i-1)*(i-1)+j;
      mC_l(ind)=mCon(j);
      mC_un(ind)=mxun(j);
      duC_l(ind)=dCon(j);
    end;
  end
end
