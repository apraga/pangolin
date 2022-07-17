%
% advection nord-sud par secteur
%
function [mC_l,mC_un,duC_l] = adv_u5(C_l,u_s,mco,dco,nlat,nlon,dt)
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
      dCon(j)=(Con(j+1)-Con(j-1))*0.5/dlon;
      dC_min=2*min(abs(Con(j+1)-Con(j)),abs(Con(j-1)-Con(j)))/dlon;
      ypm=min(abs(dCon(j)),dC_min);
      %    extr=sign(Con(j+1)-Con(j))*sign(-Con(j-1)+Con(j));
      %    extr=(extr+1)/2;
      extr=1.;
      dCon(j)=sign(dCon(j))*ypm*extr;
    end;
    dCon(1)=(Con(2)-Con(jl))*0.5/dlon;
    dC_min=2*min(abs(Con(2)-Con(1)),abs(Con(jl)-Con(1)))/dlon;
    ypm=min(abs(dCon(1)),dC_min);
    %    extr=sign(Con(2)-Con(1))*sign(-Con(jl)+Con(1));
    %    extr=(extr+1)/2;
    dCon(1)=sign(dCon(1))*ypm*extr;
    %
    dCon(jl)=(Con(1)-Con(jl-1))*0.5/dlon;
    dC_min=2*min(abs(Con(1)-Con(jl)),abs(Con(jl-1)-Con(jl)))/dlon;
    ypm=min(abs(dCon(jl)),dC_min);
    %    extr=sign(Con(1)-Con(jl))*sign(-Con(jl-1)+Con(jl));
    %    extr=(extr+1)/2;
    dCon(jl)=sign(dCon(jl))*ypm*extr;
    %
    %
    for j=2:jl;
      am=dCon(j-1);
      bm=Con(j-1)-am*0.5*dlon;
      ap=dCon(j);
      bp=Con(j)-ap*0.5*dlon;
      up=max(u_s(i,j),0.);
      cd=up*dt;
      f_u(i,j)=am*0.5*(2*dlon-cd)*cd+bm*cd;
      f_u1(i,j)=cd;
      un=min(u_s(i,j),0.);
      cd=-un*dt;
      f_u(i,j)=f_u(i,j)-ap*cd*cd*0.5-bp*cd;
      f_u1(i,j)=f_u1(i,j)-cd;
    end;
    j=1;
    am=dCon(jl);
    bm=Con(jl)-am*0.5*dlon;
    ap=dCon(j);
    bp=Con(j)-ap*0.5*dlon;
    up=max(u_s(i,j),0.);
    cd=up*dt;
    f_u(i,j)=am*0.5*(2*dlon-cd)*cd+bm*cd;
    f_u1(i,j)=cd;
    un=min(u_s(i,j),0.);
    cd=-un*dt;
    f_u(i,j)=f_u(i,j)-ap*cd*cd*0.5-bp*cd;
    f_u1(i,j)=f_u1(i,j)-cd;

    %
    for j=1:jl-1;
      dx(j)=(f_u(i,j+1)-f_u(i,j))*dlat;
      dxun(j)=(f_u1(i,j+1)-f_u1(i,j))*dlat;
      mCon(j)=-dx(j)+Con(j)*dlat*dlon;
      mxun(j)=-dxun(j)+xun(j)*dlat*dlon;
      %
    end;
    j=jl;
    dx(j)=(f_u(i,1)-f_u(i,j))*dlat;
    dxun(j)=(f_u1(i,1)-f_u1(i,j))*dlat;
    mCon(j)=-dx(j)+Con(j)*dlat*dlon;
    mxun(j)=-dxun(j)+xun(j)*dlat*dlon;
    %
    % Calcul des dérivées
    %
    for j=1:jl;
      Con(j)=mCon(j)/mxun(j);
    end;

        %
    for j=2:jl-1;
      dCon(j)=(Con(j+1)-Con(j-1))*0.5/dlon;
      %        extr=sign(Con(j+1)-Con(j))*sign(-Con(j-1)+Con(j));
      %        extr=(extr+1)/2;
      %        dCon(j)=dCon(j)*extr;
    end;
    dCon(1)=(Con(2)-Con(jl))*0.5/dlon;
    dCon(jl)=(Con(1)-Con(jl-1))*0.5/dlon;
    %
    %dCon=dCon*0.;
    %
    for j=1:jl;
      ind=3*(i-1)*(i-1)+j;
      mC_l(ind)=mCon(j);
      mC_un(ind)=mxun(j);
      duC_l(ind)=dCon(j);
    end;
  end;
