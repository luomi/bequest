clear
clc
%%

data=[34050	35850	28000	47800	29000	15000;
    175000	159050	143000	183000	130000	116200;
    515000	445000	395000	470000	356510	326000]';

T=94;
periods=T-65+1;

W=eye(periods/5*3);

Data_Mom=data(:);

Ind=20000;

r=.05;
% sigma=3;
betta=.975;
% delta=3;
% omega_bar=93.7;
% gamma=3;
% phi= 12.06;
p_0=[5.0804,0.1768,0.0853,5.0804,17.1028,52.6746,5.79,2.2,10]'; % median from SSQs
gender = 1;
% p_0=   1.0e+003 * [ 5.790000000000000   2.200000000000000                   0 ]
% C_F=parameter(1);
% LTC_pc=parameter(2);
% 
% chi_LTC=parameter(3);  
% 
% omega_bar=parameter(4);
% gamma= parameter(5);
% phi= parameter(6);
lb=[0,0,0,0,0,0,0,0,0];
ub=[1000,1000,1000,1000,1000,1000,1000,1000,1000];


%% Shocks;
rand('seed',0);
%
HS_shock=rand(periods,Ind);
rand('seed',1e30);
%
HC_shock=rand(periods,Ind);
rho=.6;
mu=[0;0];
sigma1=[1,rho;rho,1];
randn('seed',1e5);
%
x=mvnrnd(mu,sigma1,Ind);

rand('seed',5e5);
%
In_Health=rand(Ind,1);

%%
Est_F_only_parameters=g(p_0,T,gender,Data_Mom,W,r,betta,HS_shock,HC_shock,x,In_Health)
% Est_F_only_parameters=@(p) g(p,T,gender,Data_Mom,W,r,betta,HS_shock,HC_shock,x,In_Health);
% eta_0=.6; eta_1=.8; gamma_0=.7; gamma_1=1.5; delta_0=.5; 
% weights=[1,1,1,1,1,1,1,1,1]; ParallelNodes=40;

% p_out= Pounders_With_Reg(Est_F_only_parameters,p_0,1,eta_0,eta_1,gamma_0,gamma_1,delta_0,lb,ub,weights,ParallelNodes)
