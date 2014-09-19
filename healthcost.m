function [cost,prob]=healthcost(s,j, h,sex)


age=j;
iid_hc_st=h;
h_state=s;
%%
m=3;
norm_grid=[(-.441-1.96)/2,0,(.441+.96/2)];
prob_grid=ones(1,3)/m;
if h_state<4
%%
b_age=2.300422172;
b_age2=-.0498040356;
b_age3=.0004693922;
b_age4=-1.62241e-6;
b_h2=.0627911333;
b_h3=-1.026224783;
b_male=-.1218282043;
b_agemale=.0009363926;
b_heal2_age=.000825404;
b_heal3_age=.016834925;
b_cons=-39.32217509;

mean=b_cons+b_age*age+b_age2*(age^2)+b_age3*(age^3)+b_age4*(age^4)...
    +b_h2*(h_state==2)+b_h3*(h_state==3)...
    +b_male*(sex==1)+b_agemale*(sex==1)*age+...
    b_heal2_age*(h_state==2)*age+b_heal3_age*(h_state==3)*age;


v_age=2.558199073;
v_age2=-.0620090455;
v_age3=.0006314104;
v_age4=-2.29509e-6;
v_h2=1.03720271;
v_h3=2.077973606;
v_male=1.881918378;
v_agemale=-.0252754256;
v_heal2_age=-.0094832478;
v_heal3_age=-.0195496742;
v_cons=-35.03753481;

var=v_cons+v_age*age+v_age2*(age^2)+v_age3*(age^3)+v_age4*(age^4)...
    +v_h2*(h_state==2)+v_h3*(h_state==3)...
    +v_male*(sex==1)+v_agemale*(sex==1)*age...
    +v_heal2_age*(h_state==2)*age+v_heal3_age*(h_state==3)*age;
sd_est=max(var,0)^.5;

l_cost=mean+sd_est*norm_grid(iid_hc_st);
cost=exp(l_cost);


prob=prob_grid(iid_hc_st);
% cost=0;
elseif h_state==4
   cost=4;
   prob_grid=[1,0,0];
   prob=prob_grid(iid_hc_st);
end

