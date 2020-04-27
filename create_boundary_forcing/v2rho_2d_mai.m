function [var_rho]=v2rho_2d_mai(var_v)

[L,Mp]=size(var_v);
Lp=L+1;
Lm=L-1;
var_rho=zeros(Lp,Mp);
var_rho(2:L,:)=0.5*(var_v(1:Lm,:)+var_v(2:L,:));
var_rho(1,:)=var_rho(2,:);
var_rho(Lp,:)=var_rho(L,:);

return

