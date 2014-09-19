function moments=Sim_Moments(parameter,e_states,a_grid,T,r,data)

%% Define data
age_0=data(:,1);
wealth_0=data(:,2);
income_quint=data(:,3);
initial_health=data(:,4);
gender=data(:,5);
assign_vars=data(:,6:end);
[~,n_type_var]=size(assign_vars);
ma = length(a_grid);


[Ind,~]=size(data);

gender_in=1;
Final_age=55;
Period=T-Final_age+1;
%%
rng('default');
rng(0);
HS_shock=rand(Period,Ind);
% rand('seed',1e30);
rng(2^30);
HC_shock=rand(Period,Ind);
% rand('seed',5e5);
rng(5e5);
p_assign_rv=rand(Ind,1);
n_type_var=n_type_var+1;

%% Generate policy grids


tic
        
[V_sj1_1, P_sj1_1,C_sj1_1,CO_sj1_1,CF_sj1_1]=EGM_Simplified(Final_age,1,T,parameter((n_type_var+1):(n_type_var+9)),gender_in,r,a_grid,e_states);
[V_sj2_1, P_sj2_1,C_sj2_1,CO_sj2_1,CF_sj2_1]=EGM_Simplified(Final_age,2,T,parameter((n_type_var+1):(n_type_var+9)),gender_in,r,a_grid,e_states);
[V_sj3_1, P_sj3_1,C_sj3_1,CO_sj3_1,CF_sj3_1]=EGM_Simplified(Final_age,3,T,parameter((n_type_var+1):(n_type_var+9)),gender_in,r,a_grid,e_states);
[V_sj4_1, P_sj4_1,C_sj4_1,CO_sj4_1,CF_sj4_1]=EGM_Simplified(Final_age,4,T,parameter((n_type_var+1):(n_type_var+9)),gender_in,r,a_grid,e_states);
[V_sj5_1, P_sj5_1,C_sj5_1,CO_sj5_1,CF_sj5_1]=EGM_Simplified(Final_age,5,T,parameter((n_type_var+1):(n_type_var+9)),gender_in,r,a_grid,e_states);


[V_sj1_2, P_sj1_2,C_sj1_2,CO_sj1_2,CF_sj1_2]=EGM_Simplified(Final_age,1,T,parameter((n_type_var+10):(n_type_var+18)),gender_in,r,a_grid,e_states);
[V_sj2_2, P_sj2_2,C_sj2_2,CO_sj2_2,CF_sj2_2]=EGM_Simplified(Final_age,2,T,parameter((n_type_var+10):(n_type_var+18)),gender_in,r,a_grid,e_states);
[V_sj3_2, P_sj3_2,C_sj3_2,CO_sj3_2,CF_sj3_2]=EGM_Simplified(Final_age,3,T,parameter((n_type_var+10):(n_type_var+18)),gender_in,r,a_grid,e_states);
[V_sj4_2, P_sj4_2,C_sj4_2,CO_sj4_2,CF_sj4_2]=EGM_Simplified(Final_age,4,T,parameter((n_type_var+10):(n_type_var+18)),gender_in,r,a_grid,e_states);
[V_sj5_2, P_sj5_2,C_sj5_2,CO_sj5_2,CF_sj5_2]=EGM_Simplified(Final_age,5,T,parameter((n_type_var+10):(n_type_var+18)),gender_in,r,a_grid,e_states);

toc

% gender_in=0;
% [V_sj1_1_f, P_sj1_1_f,C_sj1_1_f]=EGM_Simplified(55,1,T,parameter((n_type_var+1):(n_type_var+8)),fixed_parameter,gender_in,r,betta);
% [V_sj2_1_f, P_sj2_1_f,C_sj2_1_f]=EGM_Simplified(55,2,T,parameter((n_type_var+1):(n_type_var+8)),fixed_parameter,gender_in,r,betta);
% [V_sj3_1_f, P_sj3_1_f,C_sj3_1_f]=EGM_Simplified(55,3,T,parameter((n_type_var+1):(n_type_var+8)),fixed_parameter,gender_in,r,betta);
% [V_sj4_1_f, P_sj4_1_f,C_sj4_1_f]=EGM_Simplified(55,4,T,parameter((n_type_var+1):(n_type_var+8)),fixed_parameter,gender_in,r,betta);
% [V_sj5_1_f, P_sj5_1_f,C_sj5_1_f]=EGM_Simplified(55,5,T,parameter((n_type_var+1):(n_type_var+8)),fixed_parameter,gender_in,r,betta);
% 
% 
% [V_sj1_2_f, P_sj1_2_f,C_sj1_2_f]=EGM_Simplified(55,1,T,parameter((n_type_var+9):(n_type_var+16)),fixed_parameter,gender_in,r,betta);
% [V_sj2_2_f, P_sj2_2_f,C_sj2_2_f]=EGM_Simplified(55,2,T,parameter((n_type_var+9):(n_type_var+16)),fixed_parameter,gender_in,r,betta);
% [V_sj3_2_f, P_sj3_2_f,C_sj3_2_f]=EGM_Simplified(55,3,T,parameter((n_type_var+9):(n_type_var+16)),fixed_parameter,gender_in,r,betta);
% [V_sj4_2_f, P_sj4_2_f,C_sj4_2_f]=EGM_Simplified(55,4,T,parameter((n_type_var+9):(n_type_var+16)),fixed_parameter,gender_in,r,betta);
% [V_sj5_2_f, P_sj5_2_f,C_sj5_2_f]=EGM_Simplified(55,5,T,parameter((n_type_var+9):(n_type_var+16)),fixed_parameter,gender_in,r,betta);

gender(gender==0)=1;
P_sj1_1_f=P_sj1_1;
P_sj2_1_f=P_sj2_1;
P_sj3_1_f=P_sj3_1;
P_sj4_1_f=P_sj4_1;
P_sj5_1_f=P_sj5_1;

P_sj1_2_f=P_sj1_2;
P_sj2_2_f=P_sj2_2;
P_sj3_2_f=P_sj3_2;
P_sj4_2_f=P_sj4_2;
P_sj5_2_f=P_sj5_2;

tic
%%





period=zeros(Ind,1);

%% Assign Preference types
beta_c=parameter(1);
beta_assign=parameter(2:n_type_var);
pref_type=zeros(Ind,1);

for ii=1:Ind
    type_1_prob=exp(beta_c+beta_assign'*(assign_vars(ii,:))')/(1+exp(beta_c+beta_assign'*(assign_vars(ii,:))'));
    pref_type(ii)=1+(p_assign_rv(ii)>type_1_prob);
    period(ii)=T-age_0(ii)+1;
end






%%

Inc_state=zeros(Ind,1); %Individual income level
IW_state=zeros(Ind,1);


HS_state=cell(Ind,1); %Realized individual health state
HC_state=cell(Ind,1); %Realized individual Health cost state
path=cell(Ind,1); %Sequence of health shocks
cost_path=cell(Ind,1); %Sequence of health cost shocks
need_path=cell(Ind,1); %Sequence of need state shocks
deathtime=zeros(Ind,1);

Simulated_Assets=cell(Ind,1);

assets_t=cell(Period,1);


HC_cost=zeros(4,114,3);
HC_Prob=zeros(4,114,3);
HC_Prob_f=HC_Prob;
HC_cost_f=HC_cost;

mu_e    = [parameter(n_type_var+6);parameter(n_type_var+15)];
sigma_e = [exp(parameter(n_type_var+9));exp(parameter(n_type_var+18))];

epsilon_grid1    = nodeunif(e_states,mu_e(1)-3.*sigma_e(1),mu_e(1)+3.*sigma_e(1));
epsilon_grid2    = nodeunif(e_states,mu_e(2)-3.*sigma_e(2),mu_e(2)+3.*sigma_e(2));
epsilon_prob1    = normpdf(epsilon_grid1,mu_e(1),sigma_e(1));
epsilon_prob1    = epsilon_prob1/sum(epsilon_prob1);
epsilon_prob2    = normpdf(epsilon_grid2,mu_e(2),sigma_e(2));
epsilon_prob2    = epsilon_prob2/sum(epsilon_prob2);

for ii = 1:Ind
    if pref_type(ii) == 1
        rng(1234);
        uni1     = rand(Period,1);
        cumprob1 = [0; cumsum(epsilon_prob1)];
        for ee = 1:e_states
            need_path{ii}((uni1 > cumprob1(ee)) & (uni1 <= cumprob1(ee+1)))   = ee;
        end
    elseif pref_type(ii) == 2
        rng(1234);
        uni2     = rand(Period,1);
        cumprob2 = [0; cumsum(epsilon_prob2)];
        for ee = 1:e_states
            need_path{ii}((uni2 > cumprob2(ee)) & (uni2 <= cumprob2(ee+1)))   = ee;
        end
    end
end

%% Shocks;
for ii=1:4
    for jj=54:114
        for kk=1:3
            [HC_cost(ii,jj,kk),HC_Prob(ii,jj,kk)]=healthcost(ii,jj,kk,1);
            [HC_cost_f(ii,jj,kk),HC_Prob_f(ii,jj,kk)]=healthcost(ii,jj,kk,0);
        end
    end
end

% HC_Index=55:T;
% HC_prob_cutoffs=cumsum(HC_Prob(HC_Index,:,:),3);
%%
Tran_mat=cell(Period,2);
for ii=1:Period
    Tran_mat{ii,1}=Health_Transition(55+ii-1,1);
    Tran_mat{ii,2}=Health_Transition(55+ii-1,2);
%     HS_Prob(ii,1,:)=Tran_mat{ii}(1,:);
%     HS_Prob(ii,2,:)=Tran_mat{ii}(2,:);
%     HS_Prob(ii,3,:)=Tran_mat{ii}(3,:);
%     HS_Prob(ii,4,:)=Tran_mat{ii}(4,:);
end
% HS_prob_cutoffs=cumsum(HS_Prob,3);


%%
% for i=1:Ind
% %     shock2=HS_shock(:,i);
% %     health_1=(repmat(shock2,1,4)>HS_prob_cutoffs(:,:,1));
%     health_2=(repmat(shock2,1,4)>HS_prob_cutoffs(:,:,2));
%     health_3=(repmat(shock2,1,4)>HS_prob_cutoffs(:,:,3));
%     HS_state{i}=health_1+health_2+health_3+ones(size(health_1));
%     
%     
%     shock1=HC_shock(:,i);
%     cost_1=(repmat(shock1,1,4)>HC_prob_cutoffs(:,:,1));
%     cost_2=(repmat(shock1,1,4)>HC_prob_cutoffs(:,:,2));
%     HC_state{i}=cost_1+cost_2+ones(size(cost_1));
    
%     path{i}(1)=initial_health;
%     
%     Inc_state(i)=income_quint(i);
%     IW_state(i)=wealth_0(i);
% end



for i=1:Ind
    path{i}(1)=initial_health(i);
    sex_i=1+(gender(i)==0);
    if sex_i==1
        HC_Prob_i=HC_Prob;
    elseif sex_i==2
        HC_Prob_i=HC_Prob_f;
    end
    
    Inc_state(i)=income_quint(i);
    IW_state(i)=wealth_0(i);
    
    for k=1:(period(i))
%         path{i}(k+1)=HS_state{i}(k,path{i}(k));
%         cost_path{i}(k)=HC_state{i}(k,path{i}(k));
        path{i}(k+1)=sum((cumsum(Tran_mat{age_0(i)-55+k,sex_i}(path{i}(k),:))<repmat(HS_shock(k,i),1,4)))+1;
        cost_path{i}(k)=sum((squeeze(cumsum(HC_Prob_i(path{i}(k),age_0(i)+k-1,:)))'<...
           repmat(HC_shock(k,i),1,length(HC_Prob_i(path{i}(k),age_0(i)+k-1,:)))))+1;
    end
    path{i}(period(i)+1)=[];
    
    if isempty(find(path{i}==4,1))==0
        deathtime(i)=find(path{i}==4,1);
        %Since path{i}(2)=state at time 66,path{i}(3)=state at age 67, we
        %see that path{i}=state at time 65+i-1.  Thus, the agents can at
        %most live to time 100, corresponding to path{i}(36).  Thus, a
        %deathtime(i)=k corresponds to a death at time 65+k-1.
    elseif isempty(find(path{i}==4,1))==1
        deathtime(i)=period(i);
    end
    
end


for i=1:Ind
    start_age=age_0(i)-55+1;
    Simulated_Assets{i}(1:(start_age-1))=0/0;
    Simulated_Assets{i}(start_age)=IW_state(i);  %Assets you start with at age 55
    
    if Inc_state(i)==1% Psj is the policy matrix.  The rows denote health states, while the columns denote health cost states;
        if pref_type(i)==1
            if gender(i)==1
                P_sj=P_sj1_1;
            elseif gender(i)==0
                P_sj=P_sj1_1_f;
            end
        elseif pref_type(i)==2
            if gender(i)==1
                P_sj=P_sj1_2;
            elseif gender(i)==0
                P_sj=P_sj1_2_f;
            end
        end

        
    elseif Inc_state(i)==2
        if pref_type(i)==1
            if gender(i)==1
                P_sj=P_sj2_1;
            elseif gender(i)==0
                P_sj=P_sj2_1_f;
            end
        elseif pref_type(i)==2
            if gender(i)==1
                P_sj=P_sj2_2;
            elseif gender(i)==0
                P_sj=P_sj2_2_f;
            end
        end
        
        
    elseif Inc_state(i)==3
        if pref_type(i)==1
            if gender(i)==1
                P_sj=P_sj3_1;
            elseif gender(i)==0
                P_sj=P_sj3_1_f;
            end
        elseif pref_type(i)==2
            if gender(i)==1
                P_sj=P_sj3_2;
            elseif gender(i)==0
                P_sj=P_sj3_2_f;
            end
        end
        
    elseif Inc_state(i)==4
        if pref_type(i)==1
            
            if gender(i)==1
                P_sj=P_sj4_1;
            elseif gender(i)==0
                P_sj=P_sj4_1_f;
            end
            
        elseif pref_type(i)==2
            if gender(i)==1
                P_sj=P_sj4_2;
            elseif gender(i)==0
                P_sj=P_sj4_2_f;
            end
        end;
    elseif Inc_state(i)==5
        if pref_type(i)==1
            if gender(i)==1
                P_sj=P_sj5_1;
            elseif gender(i)==0
                P_sj=P_sj5_1_f;
            end
        elseif pref_type(i)==2
            if gender(i)==1
                P_sj=P_sj5_2;
            elseif gender(i)==0
                P_sj=P_sj5_2_f;
            end
        end;
    end
    
    
    
    i_path=path{i};
    i_cost_path=cost_path{i};
    
    for k=1:(deathtime(i)-1)
        ee = need_path{i}(k);
        P_sj_temp = P_sj{i_path(k),i_cost_path(k),start_age-1+k}(((ee-1)*ma+1):(ee*ma));
        Simulated_Assets{i}(start_age+k)=interp1(a_grid, P_sj_temp, Simulated_Assets{i}(start_age+k-1),[],'extrap');
                %%Check the removal of k+1 to k in the line above in the P_sj indexing.... This
                %%comes due to the fact that I recoded the model sulution
                %%to ensure thaat the first policy function is in the first
                %%cell of the object.
        %The below discussion is for me to remind myself about indices.

        
        %This calculates the assets you hold at 65+k-1.  Thus,
        %Simulate_Assets{i}(k+1)=assets held at start of period 65+k.
        %To calculate this, we take the preious periods asset holdings,
        %Simulate_Assets{i}(65+k), or the asset holdings one held at
        %65+k-1.  Thus, when k=1, then this is simply the assets you
        %started with at time t=65.  You then calculate the assets you held
        %at time 66.
        
        %To do this, we interpolate over our the policy funcitons we
        %calculated.  Note that they way P_sj is saved, P_sj(:,:,2)
        %corresponds to the policy function at t=65, P_sj(:,:,3)
        %corresponds to the policy function at t=66, and in general
        %P_sj(:,:,m) corresponds to the policy function at t=65+m-2.
        %Therefore, to calculate Simulated_Assets{i}(65+1), we use the
        %policy function P_sj(1+1) given starting assets
        %Simulated_assets{i}(65).
        
        %The health state for time 65 is interp_path(1), 66 is interp_path(2).  Thus, for
        %Asset holdings at time 65+k (or Simulated_Assets{i}(k+1)) one
        %should interpolate over interp_path(k).  The cost path,
        %interp_cost_path(k), is exactly the same.
        
        
        
    end
    
%     bequest(i)=max((1+r)*Simulated_Assets{i}(deathtime(i)-1)-HC_cost(54+deathtime(i), i_path(deathtime(i)),i_cost_path(deathtime(i))),0);

end
for io=1:Ind
    assets_t{1}(io)=Simulated_Assets{io}(1);
end

   
for k=2:period;
        
    living_I=find(deathtime>=(start_age+k-1));
    for j=living_I'
        assets_t{k}=[assets_t{k};Simulated_Assets{j}(k)];
    end
    

end
assets_t{1}=assets_t{1}';
%% Construct Moments (FIX)
N_mom=12;
mean_vec_3yr=zeros(N_mom,1);
sd_vec_3yr=zeros(N_mom,1);
Per_25_vec_3yr=zeros(N_mom,1);
Per_50_vec_3yr=zeros(N_mom,1);
Per_75_vec_3yr=zeros(N_mom,1);

for jj=1:N_mom
    mean_vec_3yr(jj)=nanmean([assets_t{(jj-1)+1};assets_t{(jj-1)+2};assets_t{(jj-1)+3}]);
    sd_vec_3yr(jj)=nanstd([assets_t{(jj-1)+1};assets_t{(jj-1)+2};assets_t{(jj-1)+3}]);
    Per_25_vec_3yr(jj)=prctile([assets_t{(jj-1)+1};assets_t{(jj-1)+2};assets_t{(jj-1)+3}],25);
    Per_50_vec_3yr(jj)=prctile([assets_t{(jj-1)+1};assets_t{(jj-1)+2};assets_t{(jj-1)+3}],50);
    Per_75_vec_3yr(jj)=prctile([assets_t{(jj-1)+1};assets_t{(jj-1)+2};assets_t{(jj-1)+3}],75);    
end

moments=[Per_25_vec_3yr; Per_50_vec_3yr; Per_75_vec_3yr]
toc
