function value=bequest(b,theta_beq,gamma,kappa_beq)


b=b.*(b>0)+(-kappa_beq).*(b<0);
% omega_bar=1;
% gamma= 3;
% phi= 0;
value=theta_beq./(1-gamma).*(kappa_beq+b).^(1-gamma);