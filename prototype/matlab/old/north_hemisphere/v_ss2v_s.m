%
% Interpollation d'une grille nlat x nlon+1 vers une grille secteur
function [v_s] = v_ss2v_s(v_ss,nlat,n)
%
%
for i=1:nlat; 
    for j=1:4*i-1;
    v_s(i,j)=v_ss(i,j,n);
    end;
end;
%