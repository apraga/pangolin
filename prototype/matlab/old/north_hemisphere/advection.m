% Programme principal
% Advection d'un traceur sur une hemisphere
% avec maillage par secteurs
%
% Definition de la grille
nlat=90;
nlon=360;
% Pas de temps
dt=20;
bidim = false;
%
%%%%%%%%%%%%%%%
% Init sur une bande
C = zeros(nlat,nlon);
%%%%%%%%%%%%%%%
if (bidim)
  for i=1:nlat,
    for j=1:nlon;
      C(i,j)=i/nlat;
    end;
  end;
else
  %C(3,:)=1;
%  C(:,1:6)=1;
 C(:,1:110) = 0.5;
end

%
% Initialisation par liste et grille secteur
%
[mco,dco] = grille(nlat);
[mcoo,dcoo,isec,nindp,nindm]=sub_grille(nlat);
%
% Initialisation de la circulation
%if (bidim)
  [u,v] = circul (nlat,nlon);
%else
%  [u,v] = circul_analytical (nlat,nlon,1);
%end
[u_s] = reg2sec_u(u,mco,dco,nlat,nlon);
% Correction de u et v pour conservation de la masse
%if (bidim)
%[v_ss] = cor_massv(v,mcoo,dcoo,nlat,nlon);
%[u_ss,um,div,div2m] = cor_massu(u_s,v_ss,dcoo,nindm,nindp,nlat);
%else
% Pas de conservation de la masse
[v_ss] = no_cor_mass_v(v,mcoo,dcoo,nlat,nlon);
u_ss = u_s;
%end

nsec=1;
[C_s,dC_s]= reg2sec(C,mco,dco,nlat,nlon,nsec);
[C_l,xlat,xlon] = Cs2Cl(C_s,nlat,nsec);
%
for nsec=2:3;
  [C_s,dC_s]= reg2sec(C,mco,dco,nlat,nlon,nsec);
  [C_lx,xxlat,xxlon] = Cs2Cl(C_s,nlat,nsec);
  C_l=C_l+C_lx;
  xlat=xlat+xxlat;
  xlon=xlon+xxlon;
end;

%gnu_write(0,C_l,xlat,xlon,nlat);
%X = gnu_read(1000,nlat);
%C_l = X';

% Stockage concentration initiale et calcul de la masse
C_sto=C_l;
[m0] = masse(C_l,nlat)
% Integration temporelle 
kk=0;
%u_s(1,1:3)
for k=1:1700
%for k=1001:1700
  % Advection est-ouest de u
  [mC_l,mC_un,duC_l] = adv_u_alone(C_l,u_ss,mco,dco,nlat,nlon,dt);
  % Advection nord-sud de v
  C_l2=C_l*0.;
  for nsec=1:3;    
    [v_s] = v_ss2v_s(v_ss,nlat,nsec);
    [mC_s] = Cl2Cs(mC_l,nlat,nsec);
    [mC_un_s] = Cl2Cs(mC_un,nlat,nsec);
    [duC_s] = Cl2Cs(duC_l,nlat,nsec);
    [C_s,f_v] = adv_v_alone2(mC_s,duC_s,mC_un_s,v_s,mco,dco,mcoo,dcoo,nindp,nindm,nlat,nlon,nsec,dt);

    [C_lx,xxlat,xxlon] = Cs2Cl(C_s,nlat,nsec);
    C_l2=C_l2+C_lx;
  end;
  C_l=C_l2;
  fprintf('%f %f %f \n',C_l(1),C_l(2),C_l(3));

  k
  %if (mod(k,100) == 0)
  if (k == 500 || k == 1000 || k == 1700)
    gnu_write(k,C_l,xlat,xlon,nlat);
    [m0] = masse(C_l,nlat)
  end
end;
% Fin de l'evolution temporelle
% diagnostic masse
[mass1] = masse(C_l,nlat)

