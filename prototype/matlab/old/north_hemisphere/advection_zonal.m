% Programme principal
% Advection d'un traceur sur une hemisphere
% avec maillage par secteurs
%
% Definition de la grille
nlat=90;
nlon=360;
% Pas de temps
dt=20;
%
%%%%%%%%%%%%%%%
% Init sur une bande
C = zeros(nlat,nlon);
C(:,1:5)=1;
%%%%%%%%%%%%%%%

%
% Initialisation par liste et grille secteur
%
[mco,dco] = grille(nlat);
[mcoo,dcoo,isec,nindp,nindm]=sub_grille(nlat);

%name = strcat ('list_sub_grid.txt');
%fid= fopen (name,'w');
%for n=1:nlat;
%  for ic=1:4*n-1;
%    count= fprintf (fid,' %f %f %d %d \n',mcoo(n,ic),dcoo(n,ic),nindp(n,ic),nindm(n,ic));
%  end
%end
%fclose(fid);

% Initialisation de la circulation
[u,v] = circul_analytical (nlat,nlon,1);
[u_s] = reg2sec_u(u,mco,dco,nlat,nlon);
% Correction de u et v pour conservation de la masse
%[v_ss] = cor_massv(v,mcoo,dcoo,nlat,nlon);
%[u_ss,um,div,div2m] = cor_massu(u_s,v_ss,dcoo,nindm,nindp,nlat);
% Pas de conservation de la masse
u_ss = u_s;
v_ss = v;
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

gnu_write(0,C_l,xlat,xlon,nlat);

% Use previously computed values
%X = gnu_read(1000,nlat);
%C_l=X';
%
% Stockage concentration initiale et calcul de la masse
C_sto=C_l;
[m0] = masse(C_l,nlat)
% Integration temporelle 
dlat=90/nlat;
for k=1:250
  % Advection est-ouest de u
  [mC_l,mC_un,duC_l] = adv_u_alone(C_l,u_ss,mco,dco,nlat,nlon,dt);
  % Mise a jour des concentrations
  for i=1:nlat;
    jl=3*(2*i-1);
    coslat=cos((90-dlat/2-(i-1)*dlat)*pi/180);
    dlon=360*coslat/jl;
    dS = dlon/dlat;
    for j=1:jl;
      ind=3*(i-1)*(i-1)+j;
      C_l(ind) = mC_l(ind)/dS;
      C_l(ind) = mC_l(ind)/mC_un(ind);
    end;
  end;

  if (mod(k,10) == 0)
    gnu_write(k,C_l,xlat,xlon,nlat);
  end
  k
end;
% Fin de l'elvolution temporelle
% diagnostic masse
[mass1] = masse(C_l,nlat)

