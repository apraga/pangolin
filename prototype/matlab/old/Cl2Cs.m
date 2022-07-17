%C_liste vers C_secteur
function [C_s] = Cl2Cs(C_l,nlat,nsec);
%
%
n2=3*nlat*nlat;
for n=1:n2;
    i=fix(sqrt((n-0.1)/3))+1;
    j=n-3*(i-1)*(i-1);
    ns=fix((j-0.1)/(2*i-1))+1;
    j2=j-(ns-1)*(2*i-1);
%    if(ns < 1 or ns > 1) C_s(i,j2)=C_l(n);
    C(i,j2,ns)=C_l(n);
end;
for i=1:nlat;
    for j=1:2*i-1;
        C_s(i,j)=C(i,j,nsec);
    end;
end;
