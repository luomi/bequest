tic;
%%
r=.01;
a_grid=[(0:1:99)';  (100:100:5000)'];

Ind=20000;
period=100-65+1;

Income=[[0,.05,.25,.5,.75,.9, .95]',[0, 5, 10, 20, 30, 50, 100]'];
Initial_Wealth=[ [0,.25,.50,.75,.95]', [0,25,116,306,1154]'];
%%
HS_Prob=zeros(period,4,4);
Inc_state=zeros(Ind,1); IW_state=zeros(Ind,1);


HS_state=cell(Ind,1);
HC_state=cell(Ind,1);
path=cell(Ind,1);
cost_path=cell(Ind,1);
deathtime=zeros(Ind,1);

Simulated_Assets=cell(Ind,1);
bequest=zeros(Ind,1);
assets_t=cell(period,1);
mean_vec=zeros(period,1);
sd_vec=zeros(period,1);
Per_25_vec=zeros(period,1);
Per_50_vec=zeros(period,1);
Per_75_vec=zeros(period,1);

mean_delta_set=zeros(period,4);
sd_delta_set=zeros(period,4);
bequest_delta_set=zeros(Ind,4);
DT_delta_set=zeros(Ind,4);
delta_set=zeros(period,4);
Per_25_delta_set=zeros(period,4);
Per_50_delta_set=zeros(period,4);
Per_75_delta_set=zeros(period,4);


%% Shocks;
rand('seed',0);
HS_shock=rand(period,Ind);
rand('seed',1e30);
HC_shock=rand(period,Ind);


%Income/Wealth Correlation
rho=.6;
mu=[0;0];
sigma=[1,rho;rho,1];
randn('seed',1e5);
x=mvnrnd(mu,sigma,Ind);
x1=normcdf(x(:,1)); x2=normcdf(x(:,2));


IW_shock=x1;
Income_shock=x2;
rand('seed',5e5);
In_Health=rand(Ind,1);

for delta_set=1:4
    if delta_set==1;
        load('delta_dot4.mat');
    elseif delta_set==2;
        load('delta_1dot3.mat');
    elseif delta_set==3;
        load('delta_4dot.mat');
    elseif delta_set==4;
        load('delta_8dot.mat');
    end
    
    [HC_cost,HC_Prob]=healthcost(:,:,:);
    HC_Index=65:100;
    HC_prob_cutoffs=cumsum(HC_Prob(HC_Index,:,:),3);
    %%
    Tran_mat=cell(1,period);
    for i=1:period
        Tran_mat{i}=Health_Transition(65+i);
        HS_Prob(i,1,:)=Tran_mat{i}(1,:);
        HS_Prob(i,2,:)=Tran_mat{i}(2,:);
        HS_Prob(i,3,:)=Tran_mat{i}(3,:);
        HS_Prob(i,4,:)=Tran_mat{i}(4,:);
    end
    HS_prob_cutoffs=cumsum(HS_Prob,3);
      
    
    %%
    for i=1:Ind
        shock2=HS_shock(:,i);
        health_1=(repmat(shock2,1,4)>HS_prob_cutoffs(:,:,1));
        health_2=(repmat(shock2,1,4)>HS_prob_cutoffs(:,:,2));
        health_3=(repmat(shock2,1,4)>HS_prob_cutoffs(:,:,3));
        HS_state{i}=health_1+health_2+health_3+ones(size(health_1));
        
        
        shock1=HC_shock(:,i);
        cost_1=(repmat(shock1,1,4)>HC_prob_cutoffs(:,:,1));
        cost_2=(repmat(shock1,1,4)>HC_prob_cutoffs(:,:,2));
        HC_state{i}=cost_1+cost_2+ones(size(cost_1));
        
        path{i}(1)=1+1*(In_Health(i)>.75)+1*(In_Health(i)>.95);
        
        Inc_state(i)=Income(find(Income(:,1)<Income_shock(i),1,'last'),2);
        IW_state(i)=Initial_Wealth(find(Initial_Wealth(:,1)<IW_shock(i), 1, 'last' ),2);
    end
    
    
    
    for i=1:Ind
        for k=1:(period)
            path{i}(k+1)=HS_state{i}(k,path{i}(k));
            cost_path{i}(k)=HC_state{i}(k,path{i}(k));
        end
        path{i}(period+1)=[];
        
        if isempty(find(path{i}==4,1))==0
            deathtime(i)=find(path{i}==4,1);
            %Since path{i}(2)=state at time 66,path{i}(3)=state at age 67, we
            %see that path{i}=state at time 65+i-1.  Thus, the agents can at
            %most live to time 100, corresponding to path{i}(36).  Thus, a
            %deathtime(i)=k corresponds to a death at time 65+k-1.
        elseif isempty(find(path{i}==4,1))==1
            deathtime(i)=period;
        end
        
    end
    
    
    for i=1:Ind
        Simulated_Assets{i}(1)=IW_state(i);  %Assets you start with at age 65
        if Inc_state(i)==0
            P_sj=P_sj1;% Psj is the policy matrix.  The rows denote health states, while the columns denote health cost states;
        elseif Inc_state(i)==5
            P_sj=P_sj2;
        elseif Inc_state(i)==10
            P_sj=P_sj3;
        elseif Inc_state(i)==20
            P_sj=P_sj4;
        elseif Inc_state(i)==30
            P_sj=P_sj5;
        elseif Inc_state(i)==50
            P_sj=P_sj6;
        elseif Inc_state(i)==100
            P_sj=P_sj7;
        end
        
        
        
        
        i_path=path{i};
        i_cost_path=cost_path{i};
        
        for k=1:(deathtime(i)-1)
            Simulated_Assets{i}(k+1)=interp1(a_grid, P_sj{i_path(k),i_cost_path(k),k+1}, Simulated_Assets{i}(k));
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
        
        bequest(i)=max((1+r)*Simulated_Assets{i}(deathtime(i)-1)-HC_cost(64+deathtime(i), i_path(deathtime(i)),i_cost_path(deathtime(i))),0);
        
    end
    
    for k=1:period;
        
        living_I=find(deathtime>k);
        for j=living_I'
            assets_t{k}=[assets_t{k};Simulated_Assets{j}(k)];
        end
        
        mean_vec(k)=mean(assets_t{k});
        sd_vec(k)=std(assets_t{k});
        Per_25_vec(k)=prctile(assets_t{k},25);
        Per_50_vec(k)=prctile(assets_t{k},50);
        Per_75_vec(k)=prctile(assets_t{k},75);
    end
    
    mean_delta_set(:,delta_set)=mean_vec;
    sd_delta_set(:,delta_set)=sd_vec;
    bequest_delta_set(:,delta_set)=bequest;
    DT_delta_set(:,delta_set)=deathtime;
    Per_25_delta_set(:,delta_set)=Per_25_vec;
    Per_50_delta_set(:,delta_set)=Per_50_vec;
    Per_75_delta_set(:,delta_set)=Per_75_vec;
end



toc


