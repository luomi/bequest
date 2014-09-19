function out=Utility(co, cf, alpha, sigma, delta_LTC, kappa, psi_fam)
% alpha denotes which states are active... for example, if only the long
% term care state is active, then alpha=[0,1,0]
%If  s=1 or 2 and m=1, then alpha=[1 0 1];
%
[m1,~]     = size(co);
[m2,~]     = size(cf);
% m           = max(m1,m2);
% v=0;
co2         = max(co,0);
co(isnan(co)==0)    = co2(isnan(co)==0);

cf2         = max(cf,0);
cf(isnan(cf)==0)    = cf2(isnan(cf)==0);

% Utility_Vec=zeros(m,1);
% Utility_Vec(:,1)=co.^(1-sigma)./(1-sigma);
% Utility_Vec(:,2)=delta_LTC.*(cf+kappa_LTC).^(1-sigma)./(1-sigma)-v*ones(m,1);

out         = (1+alpha.*delta_LTC).*((co.^psi_fam).*((cf+kappa).^(1-psi_fam)))...
                .^(1-sigma)./(1-sigma);
% out=(alpha*Utility_Vec')';
