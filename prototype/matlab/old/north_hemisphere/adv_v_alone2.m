%
% advection nord-sud par secteur
%
function [C_s,f_v,mC_un_s] = adv_v_alone2(mC_s,duC_s,mC_un_s,v_s,mco,dco,mcoo,dcoo,nindp,nindm,nlat,nlon,nsec,dt)
  function C_s = update_C(nlat, dlat, mC_s)
    for i=1:nlat;
      coslat=cos((90-dlat/2-(i-1)*dlat)*pi/180);
      jl=3*(2*i-1);    
      dlon=360*coslat/jl;
      dS = dlat*dlon;

      for j=1:2*i-1;
        C_s(i,j)=mC_s(i,j)/dS;
      end;
    end;
  end

  dlat=90/nlat;
  pi=2*acos(0.);
  %
  %C_s = update_C(nlat,dlat,mC_s);
  for i=1:nlat;
    for j=1:2*i-1;
      C_s(i,j)=mC_s(i,j)/mC_un_s(i,j);
    end;
  end;


  %  Limitation de la pente 
  for i=2:nlat-1;
    for j=1:2*i-1;
      coo1=fix((mco(i,j))/dco(i-1));
      dcoo1=mco(i,j)-coo1*dco(i-1);
      coslat=cos((90-dlat/2-(i-2)*dlat)*pi/180);
      jl=2*(i-1)-1;
      dlon=120*coslat/jl;
      dcoo1=dcoo1*dlon;
      coo1=coo1+1;
      C_sm=C_s(i-1,coo1)+(dcoo1-dlon/2)*duC_s(i-1,coo1);

      %
      %     
      coo2=fix((mco(i,j))/dco(i+1));
      dcoo2=mco(i,j)-coo2*dco(i+1);
      coslat=cos((90-dlat/2-(i)*dlat)*pi/180);
      jl=2*(i+1)-1;
      dlon=120*coslat/jl;
      dcoo2=dcoo2*dlon;
      coo2=coo2+1;
      C_sp=C_s(i+1,coo2)+(dcoo2-dlon/2)*duC_s(i+1,coo2);

      %  Limitation de la pente 
      dC_min=2*min(abs(C_sp-C_s(i,j)),abs(C_sm-C_s(i,j)))/dlat;
      %
      dC_s(i,j)=(C_sp-C_sm)*0.5/dlat;
      ypm=min(abs(dC_s(i,j)),dC_min);
      dC_s(i,j)=sign(dC_s(i,j))*ypm;
      % % Lax-Wendroff
      %   dC_s(i,j)=0.5*(C_sm - C_sp)/dlat;
    end;
  end;

  for j=1:2*nlat-1;
    dC_s(nlat,j)=0.;
    %
  end;
  dC_s(1,1)=0.;

  % Pour chaque noeud, on calcule l'advection entre jm (lat precedente) et 
  % jp (lat suivante)
  % Ex: le noeud a 1-5
  for i=1:nlat-1;
    coslat=cos((90-i*dlat)*pi/180);
    for j=1:4*i-1;
      % Flux de (i,jp) vers (i+1,jm)
      jm=nindm(i,j);
      am=dC_s(i+1,jm);
      bm=C_s(i+1,jm)-am*dlat*0.5;
      jp=nindp(i,j);
      ap=dC_s(i,jp);
      bp=C_s(i,jp)-ap*dlat*0.5;
      vp=max(v_s(i,j),0.);
      %
      cd=vp*dt;
      f_v(i,j)=ap*0.5*(2*dlat-cd)*cd+bp*cd;
      f_v1(i,j)=cd;
      vn=min(v_s(i,j),0.);
 
      %
      cd=-vn*dt;
      f_v(i,j)=(f_v(i,j)-am*0.5*cd*cd-bm*cd)*dcoo(i,j)*120*coslat;
      f_v1(i,j)=(f_v1(i,j)-cd)*dcoo(i,j)*120*coslat;
    end;
  end;
  for j=1:4*nlat-1;
    f_v(nlat,j)=0.;
    f_v1(nlat,j)=0.;
  end;

  for i=1:nlat-1;
    for j=1:4*i-1;
      jm=nindm(i,j);
      jp=nindp(i,j);
      % Advection vers l'equateur
      mC_s(i+1,jm)=mC_s(i+1,jm)+f_v(i,j);
      mC_s(i,jp)=mC_s(i,jp)-f_v(i,j);

      mC_un_s(i+1,jm)=mC_un_s(i+1,jm)+f_v1(i,j);
      mC_un_s(i,jp)=mC_un_s(i,jp)-f_v1(i,j);

    end;
  end;

  for i=1:nlat;
    for j=1:2*i-1;
      C_s(i,j)=mC_s(i,j)/mC_un_s(i,j);
    end;
  end;

  %C_s = update_C(nlat,dlat,mC_s);
end
