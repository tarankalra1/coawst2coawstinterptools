function var_u=rho2u_2d_mai(var_rho);

[Lp,Mp]=size(var_rho);
M=Mp-1;
var_u=0.5*(var_rho(:,1:M)+var_rho(:,2:Mp));

return