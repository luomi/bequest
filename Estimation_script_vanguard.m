clear
clc
%%

data=csvread('core_data_in.csv');
spoid_final = data(:,1);

data_trans    = csvread('core_trans_in.csv');
age     = data_trans(:,1);
trans   = data_trans(:,2)/1000;

age_mat = reshape(55:90,[3,12])';

age_mat_before=[zeros(1,3);age_mat(1:11,:)];
age_mat_next=[age_mat(2:12,:);zeros(1,3)];
age_mat=[age_mat_before,age_mat,age_mat_next];



trans_mom_25=zeros(12,1);
trans_mom_50=zeros(12,1);
trans_mom_75=zeros(12,1);

for ii=1:12
    trans_temp=trans(ismember(age, age_mat(ii,:)));
    cell_count=sum(ismember(age, age_mat(ii,:)));

 
%     trans_mom_10(ii)=prctile(trans_temp,10);
    trans_mom_25(ii)=prctile(trans_temp,25);
    trans_mom_50(ii)=prctile(trans_temp,50);
    trans_mom_75(ii)=prctile(trans_temp,75);
%     trans_mom_90(ii)=prctile(trans_temp,90);
end
Data_trans_mom=[trans_mom_25; trans_mom_50; trans_mom_75];




% data2=xlsread('P:\nyu\SSQ_Ans.csv');

% spoid2=data2(:,12);
% spoid_final=intersect(spoid1,spoid2);
[~,i_unique,~]=unique(spoid_final);

% data(ismember(spoid1,spoid_final)==0,:)=[];
% data2(ismember(spoid2,spoid_final)==0,:)=[];
% spoid1(ismember(spoid1,spoid_final)==0)=[];
% spoid2(ismember(spoid2,spoid_final)==0)=[];
% [~,I_1]=sort(spoid1);
% [~,I_2]=sort(spoid2);
[~,I_1]=sort(spoid_final);

data=data(I_1,:);
% data2=data2(I_2,:);


data(:,5)=data(:,5)/1e3;
data(data(:,5)>10000,5)=10000;
data(data(:,3)>108,3)=108;
age=data(:,3);
% data(data(:,5)>3,:)=1;
income_group=data(:,6);
assets=data(:,5);
gender=(data(:,2)==1);
health=data(:,7);
spoid=data(:,1);

Ind = 50000;



%%
% data=xlsread('P:\nyu\Asset_Profile_Formation\state_variables.csv');
% spoid=data(:,6);
% 
% 
% 
% % data2=xlsread('P:\nyu\SSQ_Ans.csv');
% % SSQ=data2(:,1:9);[~,i_unique,~]=unique(spoid);
% data=data(i_unique,:);
% data(:,3)=data(:,3)/1e3;
% data(data(:,3)>10000,:)=[];
% age=data(:,1);
% income_group=data(:,2)+1;
% assets=data(:,3);
% gender=data(:,4);
% health=data(:,5);
% spoid=data(:,6);
% spoid2=data2(:,10);
%  opts        = statset('Display','final','MaxIter',100000);
%     [groupvec_temp,centervec]    = kmeans(SSQ,2,'Replicates',1000,...      
%       'Options',opts,'distance', 'correlation');
%% Fill in groupvec: This doesn't make much sense, and is only a fill-in until we get a full sample.
% if sum(groupvec_temp==1)<sum(groupvec_temp==2)
%     groupvec_temp=-(groupvec_temp-2)+1;   
% end
% index_s1_s2=ismember(spoid,spoid2);
% index_s2_s1=ismember(spoid2,spoid);
% mdl=LinearModel.fit(data(index_s1_s2,1:5), groupvec_temp(index_s2_s1));
% groupvec=predict(mdl,data(:,1:5));
% groupvec=round(groupvec+.1);
% % groupvec(index_s1_s2)=groupvec_temp(index_s2_s1);
% 
% groupvec(isnan(groupvec))=1+(rand(sum(isnan(groupvec)),1)>.5);
  


%%

age_mat= reshape(55:90,[3,12])';

age_mat_before=[zeros(1,3);age_mat(1:11,:)];
age_mat_next=[age_mat(2:12,:);zeros(1,3)];
age_mat=[age_mat_before,age_mat,age_mat_next];



wealth_mom_25=zeros(12,1);
wealth_mom_50=zeros(12,1);
wealth_mom_75=zeros(12,1);
for ii=1:12
    asset_temp=assets(ismember(age, age_mat(ii,:)));
    cell_count=sum(ismember(age, age_mat(ii,:)));

 
%     wealth_mom_10(ii)=prctile(asset_temp,10);
    wealth_mom_25(ii)=prctile(asset_temp,25);
    wealth_mom_50(ii)=prctile(asset_temp,50);
    wealth_mom_75(ii)=prctile(asset_temp,75);
%     wealth_mom_90(ii)=prctile(asset_temp,90);
end
Data_mom=[wealth_mom_25; wealth_mom_50; wealth_mom_75];
% plot(wealth_mom_25); hold on
% plot(wealth_mom_50);
% plot(wealth_mom_75);
% plot(wealth_mom_10);
% plot(wealth_mom_90);




T=108;
periods=T-55+1;

Data_mom_temp   = Data_mom;
Data_mom_temp(Data_mom_temp==0)=1;
W=diag(1./Data_mom_temp);



% Ind=10000;

r=.05;

p_0 = [ 3.6760;
    0.0157;
         0;

    1.2115;
    1.1110;
    1.0182;
   -0.0624;
    5.5145;
   -3.0250;
    0.6793;
    3.7243;
    1.9489;
    
    1.1031;
    0.2539;
    1.4089;
   -0.1891;
   28.2049;
   36.7143;
    1.1922;
    3.1671;
    1.7465];

% p_0 = [2.4141;
%    -0.0227;
%          0;
%          
%    -1.8930;
%     1.9695;
%     0.0244;
%     0.1719;
%    38.7109;
%    55.5078;
%     0.7470;
%     2.7441;
%     0.7286;
%     
%    -0.5113;
%     1.9208;
%     0.1868;
%    -0.8281;
%    27.3828;
%    17.2266;
%     1.4908;
%     2.8574;
%     1.9878];

% p_0=[2.9069;
%     0.055;
%     0;
%     
%    -0.2883;
%     2.0608;
%     0.1177;
%     0.0324;
%   -29.9158;
%    26.4926;
%     1.2060;
%     2.9113;
%     2.3408;
%    
%    -2.2073;
%     2.2897;
%     1.7995;
%     5.5791;
%    35.8025;
%    21.3811;
%     1.7276;
%     2.4937;
%     1.5447];


% p_0(3)=0;
% order=@(val_x) floor(log(abs(val_x))./log(10));
% magnitude_p_0=order(p_0);
% magnitude_p_0(3) =0;




% 

% magnitude_p_0=order(p_0);
% p_0=p_0./(10.^magnitude_p_0);
% p_0=p_0.*(10.^magnitude_p_0);
fixed_parameter=[log(40)]; %chi_LTC



%% Shocks;
% rng('default');
% rng(0);
% HS_shock=rand(periods,Ind);
% rng(2^30);
% HC_shock=rand(periods,Ind);

% rng(5e5);
% param_assign_rv=rand(Ind,1);

rng(5e6);
sample_vals=unidrnd(length(spoid),[Ind,1]);

data_in=[age, assets, income_group, health, gender, assets, gender];
data_in=data_in(sample_vals,:);
%%
% a_grid=[(0:1:99)';  (100:100:4900)'; (5000:1000:20000)'];
a_grid=[(0:3:99)';  (100:20:500)'; (550:50:950)'; (1000:200:4800)'; (5000:1000:20000)'];
% Est_F_only_parameters=g(p_0,fixed_parameter,a_grid, T,r,HS_shock,HC_shock, ...
%     param_assign_rv, data_in , Data_mom, W)
e_states=5;
Est_F_only_parameters=@(p) g_par(p,e_states,a_grid, T,r, data_in , Data_mom, W);


eta_0=.8; eta_1=.9; gamma_0=.25; gamma_1=1.5; delta_0=.5;

% ub=[15;15;log(20-1);log(100); log(100); 100; 100; log(200); log(100); 8;  
%      log(20-1);log(100); log(100); 100; 100; log(200); log(100); 8]';
% lb =[-15;-15;log(2-1);log(1e-3); log(1e-3); -100; 0; log(1e-1); log(1e-1); 0;  
%      log(2-1);log(1e-3); log(1e-3); -100; 0; log(1e-1); log(1e-1); 0]';  
%  
ub=[4; 1e-1; 0;
     log(5-1);log(10); log(8); 2; 80; 80; log(7); 4; log(20);
     log(5-1);log(10); log(8); 2; 80; 80; log(7); 4; log(20)];
lb=[2; -1e-1; 0;
     log(1.5-1);log(1); log(1); -2; -10; -10; log(2); 3; log(1);  
     log(1.5-1);log(1); log(1); -2; -10; -10; log(2); 3; log(1)];
 
 
% ub=ub'./(10.^magnitude_p_0);
% lb=lb'./(10.^magnitude_p_0);
% weights=[1,1 ...
%     ,1,.8,.8,.5,.5,.1,.1,.2,...
%      1,.8,.8,.5,.5,.1,.1,.2]; 

%%
% test_vec=rand(50,18).*repmat((ub-lb),50,1)+repmat(lb,50,1);
% test_vec(1,:)=[-10;10;...
%     log(4-1);log(1); log(1); -10; 20; log(5); log(5); 3.5;  
%      log(4-1);log(1); log(4); -10; 20; log(5); log(5); 3.5];  
%  
%  matlabpool(51)
%  moment_vec=cell(50,1);
% parfor ii=4:50
%     ii
%     moment_vec{ii}=Sim_Moments(test_vec(ii,:),fixed_parameter,a_grid,T,r,data_in )
% end
% weights=[1,1 ,1,.8,.8,.5,.5,.1,.1,2, 1,.8,.8,.5,.5,.1,.1,2]; 

%%

p_0(3)=0;
% ParallelNodes=51;
% matlabpool(ParallelNodes)
matlabpool close force local
matlabpool('open','local',2);
% options=psoptimset('UseParallel', 'always', 'Vectorized', 'Off', 'CompletePoll', 'Off', ...
%     'Display', 'iter', 'PollingOrder','Random', 'PollMethod', 'GSSPositiveBasis2N',...
%     'MeshAccelerator', 'On', 'Scale', 'On', 'Cache', 'On', 'CacheTol', 1e-2,...
%     'TolFun', 1e-4, 'TolX', 1e-4);
% 
%     
% % parameter_est=simulannealbnd(Est_F_only_parameters, p_0,lb,ub,options);
%  parameter_est=   patternsearch(Est_F_only_parameters, p_0,[],[],[],[],lb,ub,[],options);
%     parameter_est_final=parameter_est.*(10.^magnitude_p_0);
% % p_out= Pounders_With_Reg(Est_F_only_parameters,p_0,1,eta_0,eta_1,gamma_0,gamma_1,delta_0,lb,ub,weights,ParallelNodes)

options     = optimset('Display','iter','MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-3,'TolX',1e-3);
% problem     = createOptimProblem('fminunc','objective',Est_F_only_parameters,...
%                                  'x0',p_0,'options',options);
% ms          = MultiStart('UseParallel','always','Display','iter','TolFun',1e-3,'TolX',1e-3,'StartPointsToRun','bounds');
% [parameter_est,fval,exitflag,output,solutions] = run(ms,problem,12)
[parameter_est,fval,exitflag,output] = fminsearch(Est_F_only_parameters,p_0,options);
% parameter_est_final     = parameter_est.*(10.^magnitude_p_0)

matlabpool close
save Struct_est_PS_8


%% Plots
agevec  = linspace(56,89,12)';
plot(agevec,assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(25:36),'b+:','LineWidth',2);hold on;
legend('Model','Data','Location','Northwest');
plot(agevec,assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(13:24),'b+:','LineWidth',2);hold on;
plot(agevec,assets_moments(1:12),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(1:12),'b+:','LineWidth',2);xlabel('Age');ylabel('Assets in thousands');hold off;

agevec  = linspace(56,89,12)';
% plot(agevec,assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(25:36),'b+:','LineWidth',2);hold on;
legend('Data','Location','Northwest');
% plot(agevec,assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(13:24),'b+:','LineWidth',2);hold on;
% plot(agevec,assets_moments(1:12),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(1:12),'b+:','LineWidth',2);
xlabel('Age');ylabel('Assets in thousands');hold off;


agevec  = linspace(56,89,12)';
plot(agevec,cf_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,Data_trans_mom(25:36),'b+:','LineWidth',2);hold on;
legend('Model','Data','Location','Northwest');
plot(agevec,cf_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,Data_trans_mom(13:24),'b+:','LineWidth',2);hold on;
plot(agevec,cf_moments(1:12),'ro-','LineWidth',2);hold on;
plot(agevec,Data_trans_mom(1:12),'b+:','LineWidth',2);
xlabel('Age');ylabel('Transfer in thousands');
hold off;


[countsw binCentersw] = hist(assets,256);
[countst binCenterst] = hist(trans,256);
subplot(1,2,1);
plot(binCentersw,countsw,'ro-','LineWidth',2);legend('Wealth','Location','Northeast');
subplot(1,2,2);
plot(binCenterst,countst,'b+:','LineWidth',2);legend('Transfers','Location','Northeast');




assets1 = load('results_thetabeq_01.mat','assets_moments');
assets2 = load('results_thetabeq_05.mat','assets_moments');
assets3 = load('results_thetabeq_20.mat','assets_moments');
assets4 = load('results_thetabeq_40.mat','assets_moments');

agevec  = linspace(56,89,12)';
ax(1)   = subplot(2,2,1);
plot(agevec,assets1.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\theta\times 0.1$$','interpreter','latex');
hold off;

ax(2)   = subplot(2,2,2);
plot(agevec,assets2.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\theta\times 0.5$$','interpreter','latex');
hold off;

ax(3)   = subplot(2,2,3);
plot(agevec,assets3.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\theta\times 2$$','interpreter','latex');
hold off;

ax(4)   = subplot(2,2,4);
plot(agevec,assets4.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\theta\times 4$$','interpreter','latex');
linkaxes([ax(4) ax(3) ax(2) ax(1)],'xy');
hold off;


assets1 = load('results_kbeq_01.mat','assets_moments');
assets2 = load('results_kbeq_05.mat','assets_moments');
assets3 = load('results_kbeq_20.mat','assets_moments');
assets4 = load('results_kbeq_40.mat','assets_moments');

agevec  = linspace(56,89,12)';
ax(1)   = subplot(2,2,1);
plot(agevec,assets1.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$k\times 0.1$$','interpreter','latex');
hold off;

ax(2)   = subplot(2,2,2);
plot(agevec,assets2.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$k\times 0.5$$','interpreter','latex');
hold off;

ax(3)   = subplot(2,2,3);
plot(agevec,assets3.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$k\times 2$$','interpreter','latex');
hold off;

ax(4)   = subplot(2,2,4);
plot(agevec,assets4.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$k\times 4$$','interpreter','latex');
linkaxes([ax(1) ax(2) ax(3) ax(4)],'xy');
hold off;



assets1 = load('results_mufam_01.mat','assets_moments');
assets2 = load('results_mufam_05.mat','assets_moments');
assets3 = load('results_mufam_20.mat','assets_moments');
assets4 = load('results_mufam_40.mat','assets_moments');

agevec  = linspace(56,89,12)';
ax(1)   = subplot(2,2,1);
plot(agevec,assets1.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\mu_\kappa\times 0.1$$','interpreter','latex');
hold off;

ax(2)   = subplot(2,2,2);
plot(agevec,assets2.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\mu_\kappa\times 0.5$$','interpreter','latex');
hold off;

ax(3)   = subplot(2,2,3);
plot(agevec,assets3.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\mu_\kappa\times 2$$','interpreter','latex');
hold off;

ax(4)   = subplot(2,2,4);
plot(agevec,assets4.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\mu_\kappa\times 4$$','interpreter','latex');
linkaxes([ax(3) ax(4) ax(2) ax(1)],'xy');
hold off;



assets1 = load('results_sigmafam_01.mat','assets_moments');
assets2 = load('results_sigmafam_05.mat','assets_moments');
assets3 = load('results_sigmafam_20.mat','assets_moments');
assets4 = load('results_sigmafam_40.mat','assets_moments');

agevec  = linspace(56,89,12)';
ax(1)   = subplot(2,2,1);
plot(agevec,assets1.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets1.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\sigma_\kappa\times 0.1$$','interpreter','latex');
hold off;

ax(2)   = subplot(2,2,2);
plot(agevec,assets2.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets2.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\sigma_\kappa\times 0.5$$','interpreter','latex');
hold off;

ax(3)   = subplot(2,2,3);
plot(agevec,assets3.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets3.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\sigma_\kappa\times 2$$','interpreter','latex');
hold off;

ax(4)   = subplot(2,2,4);
plot(agevec,assets4.assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,assets4.assets_moments(1:12),'ro-','LineWidth',2);grid on;
xlabel('Age');ylabel('Assets in thousands');title('$$\sigma_\kappa\times 4$$','interpreter','latex');
linkaxes([ax(3) ax(2) ax(1) ax(4)],'xy');
hold off;

