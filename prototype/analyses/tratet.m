% traces sur iso teta
LEG='SURFACE ISENTROPE';
figure (1), [C,h] = contourf (lon,lat,PTET'/100);
clabel (C,h);
xlabel ('LONGITUDE','FontSize',14);
ylabel ('LATITUDE','FontSize',14);
title({'PRESSION (mPa)',LEG},'FontSize',14);
%
figure (2), [C,h] = contourf (lon,lat,TTET');
clabel (C,h);
xlabel ('LONGITUDE','FontSize',14);
ylabel ('LATITUDE','FontSize',14);
title({'TEMPERATURE (K)',LEG},'FontSize',14);
%
figure (3), [C,h] = contourf (lon,lat,VTET');
clabel (C,h);
xlabel ('LONGITUDE','FontSize',14);
ylabel ('LATITUDE','FontSize',14);
title({'VITESSE MERIDIENNE (m/s)',LEG},'FontSize',14);
%
figure (4), [C,h] = contourf (lon,lat,UTET');
clabel (C,h);
xlabel ('LONGITUDE','FontSize',14);
ylabel ('LATITUDE','FontSize',14);
title({'VITESSE ZONALE (m/s)',LEG},'FontSize',14);
%
figure (5);
hold on;
[C,h] = contourf (lon,lat,ETET','y');
hold on
quiver(lon,lat,UTET',VTET','r');
clabel (C,h);
xlabel ('LONGITUDE','FontSize',14);
ylabel ('LATITUDE','FontSize',14);
title({'VECTEUR VITESSE (m/s) - ENERGIE CINETIQUE (m ^2 s^-^2)',LEG},'FontSize',14);
