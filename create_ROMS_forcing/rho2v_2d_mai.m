function var_v=rho2v_2d_mai(var_rho);

[Lp,Mp]=size(var_rho);
L=Lp-1;
var_v=0.5*(var_rho(1:L,:)+var_rho(2:Lp,:));

return