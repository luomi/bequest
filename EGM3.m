function [V_sj, P_sj,C_sj]=EGM3(j,y,T,parameter,r,sigma,betta_in,delta,omega_bar,gamma,phi)
% This function takes asset values a, health state s, age j, and income y,
% and returns the Value function and Policy function (asset holdings for
% period j+1) in these states.

%% 1 - Initialization of parameters
% clc
hc_states=3; % I need to make this endogenous, but for now this will work as long as it matches the healthcost.m function.
p=[ sigma %sigma
    delta %delta %   .4  1.3  4   8
    betta_in%beta
    ];
sigma=p(1);
delta=p(2);
betta=p(3);

% r=parameter(4);
C_F=parameter(1);
LTC_pc=parameter(2);

chi_LTC=parameter(3);
% 
% omega_bar=parameter(4);
% gamma= parameter(5);
% phi= parameter(6);


Final_age=j;
% final_a=a;
% final_h=h;
% final_s=s;
% 
%% Define Grids

a_grid=[(0:1:99)';  (100:100:5000)'];
m=length(a_grid);
%%
% Grids to be reported
V_sj=cell(4,hc_states,T-j+1);
P_sj=cell(3,hc_states,T-j+1);
C_sj=cell(3,hc_states,T-j+1);

% Here I am simply initializing all data structures I will need in this
% problem.  Preallocating here saves time.
I1=cell(1,hc_states);I2=cell(1,hc_states);I3=cell(1,hc_states);
Endog_1=cell(1,hc_states); Endog_2=cell(1,hc_states); Endog_3=cell(1,hc_states);
Policy_store1=cell(1,hc_states);Policy_store2=cell(1,hc_states);Policy_store3=cell(1,hc_states);
I1_skip=cell(1,hc_states);I2_skip=cell(1,hc_states);I3_skip=cell(1,hc_states);
Miss_1=cell(1,hc_states); Miss_2=cell(1,hc_states); Miss_3=cell(1,hc_states);
lower_GS1=cell(1,hc_states); lower_GS2=cell(1,hc_states);lower_GS3=cell(1,hc_states);
upper_GS1=cell(1,hc_states); upper_GS2=cell(1,hc_states);upper_GS3=cell(1,hc_states);
lower_GS1A=cell(1,hc_states); lower_GS2A=cell(1,hc_states);lower_GS3A=cell(1,hc_states);
upper_GS1A=cell(1,hc_states);  upper_GS2A=cell(1,hc_states);upper_GS3A=cell(1,hc_states);
GS_I1=cell(1,hc_states); GS_I2=cell(1,hc_states); GS_I3=cell(1,hc_states);
budget_bound=cell(1,hc_states);
Out_of_budget1=cell(1,hc_states); Out_of_budget2=cell(1,hc_states); Out_of_budget3=cell(1,hc_states);
Consumption_store1=cell(1,hc_states); Consumption_store2=cell(1,hc_states); Consumption_store3=cell(1,hc_states);
Value_store1=cell(1,hc_states); Value_store2=cell(1,hc_states); Value_store3=cell(1,hc_states);Value_store4=cell(1,hc_states);
Max_I1_T=cell(1,hc_states); Max_I2_T=cell(1,hc_states); Max_I3_T=cell(1,hc_states);
Con_fl3=cell(1,hc_states);
Gov_I1=cell(1,hc_states); Gov_I2=cell(1,hc_states); Gov_I3=cell(1,hc_states);
Con_bind3=cell(1,hc_states);
Gov_bind_I1=cell(1,hc_states); Gov_bind_I2=cell(1,hc_states); Gov_bind_I3=cell(1,hc_states);
Max_I1=cell(1,hc_states); Max_I2=cell(1,hc_states); Max_I3=cell(1,hc_states);
EYE=eye(4);
hc1=zeros(hc_states,1);hc2=hc1; hc3=hc1; hc4=hc1;
hc1_minus=zeros(hc_states,1);hc2_minus=hc1_minus; hc3_minus=hc1_minus; hc4_minus=hc1_minus;
prob1=zeros(hc_states,1);prob2=prob1; prob3=prob1; prob4=prob1;
prob1_minus=zeros(hc_states,1);prob2_minus=prob1_minus; prob3_minus=prob1_minus; prob4_minus=prob1_minus;
Value_mat1=zeros(m,hc_states); Policy_mat1=zeros(m,hc_states); Consump_mat1=zeros(m,hc_states);
Value_mat2=zeros(m,hc_states); Policy_mat2=zeros(m,hc_states); Consump_mat2=zeros(m,hc_states);
Value_mat3=zeros(m,hc_states); Policy_mat3=zeros(m,hc_states); Consump_mat3=zeros(m,hc_states);
Value_mat4=zeros(m,hc_states);
deriv_1_next=zeros(m,hc_states);
deriv_2_next=zeros(m,hc_states);
deriv_3_next=zeros(m,hc_states);
deriv_4_next=zeros(m,hc_states);

concave_L=cell(1,3);concave_H=cell(1,3);
G_NC=cell(1,3); G_C_L=cell(3,1); G_C_H=cell(1,3);
store_NC=cell(1,3); max_NC=cell(1,3);
NC_Ind=cell(1,3);
a_NC=cell(1,3);
All_Concave_I=cell(1,3);
%% T-1 Problem
% % 3a)
%

%% T-1 Problem
% % 3a)
%

for i=1:hc_states
    [hc1(i),prob1(i)]=healthcost(1,T,i);  % Here I generate the health costs and probabilities of each cost for time T.
    [hc2(i),prob2(i)]=healthcost(2,T,i);
    [hc3(i),prob3(i)]=healthcost(3,T,i);
    [hc4(i),prob4(i)]=healthcost(4,T,i);
    
    
    [hc1_minus(i),prob1_minus(i)]=healthcost(1,T-1,i); % Here I generate the health costs and probabilities of each cost for time T-1.
    [hc2_minus(i),prob2_minus(i)]=healthcost(2,T-1,i);
    [hc3_minus(i),prob3_minus(i)]=healthcost(3,T-1,i);
    [hc4_minus(i),prob4_minus(i)]=healthcost(4,T-1,i);
    
    
end
%% Here I define the value function for each state, 1,2,3.  It is the
% utility one gets from consumption, plus
% the expected value given you die next period.  We take expectation over Health
% transition  (see multiplication by Health_Transition(T-1) matrix) and
% over iid health costs for time T (see multiplication by prob1,prob2,
% etc.).


Value_Function1=@(c, a, h) Utility(c, 1, 1,[1 0 0],p)+(betta*EYE(1,:)*Health_Transition(T-1)*[prob1'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc1),0),omega_bar,gamma,phi);...
    prob2'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc2),0),omega_bar,gamma,phi);...
    prob3'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc3),0),omega_bar,gamma,phi);...
    prob4'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc4),0),omega_bar,gamma,phi)])';

Value_Function2=@(c, a, h) Utility(c, 1, 1,[1 0 0],p)+(betta*EYE(2,:)*Health_Transition(T-1)*[prob1'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc1),0),omega_bar,gamma,phi);...
    prob2'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc2),0),omega_bar,gamma,phi);...
    prob3'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc3),0),omega_bar,gamma,phi);...
    prob4'*bequest(max(repmat((1+r)*((1+r)*a'+y-c'-h),hc_states,1)-kron(ones(1,length(a)),hc4),0),omega_bar,gamma,phi)])';
%
%
Value_Function3=@(e, a, h) Utility(1, e, 1,[0 1 0],p)+(betta*EYE(3,:)*Health_Transition(T-1)*[prob1'*bequest(max(repmat((1+r)*((1+r)*a'+y-e'-h),hc_states,1)-kron(ones(1,length(a)),hc1),0),omega_bar,gamma,phi);...
    prob2'*bequest(max(repmat((1+r)*((1+r)*a'+y-e'-h),hc_states,1)-kron(ones(1,length(a)),hc2),0),omega_bar,gamma,phi);...
    prob3'*bequest(max(repmat((1+r)*((1+r)*a'+y-e'-h),hc_states,1)-kron(ones(1,length(a)),hc3),0),omega_bar,gamma,phi);...
    prob4'*bequest(max(repmat((1+r)*((1+r)*a'+y-e'-h),hc_states,1)-kron(ones(1,length(a)),hc4),0),omega_bar,gamma,phi)])';


%% T-1 Problem with Endogenous Grid
%3b) Calculate derivative of T-1 Continuation Value over asset grid.  It is
%the derivative of the bequest function if the government provided bequest
%floor is not binding, and 0 if the government bequest floor is binding.
deriv_bequest_vec=[max(((((1+r).*repmat(a_grid,1,hc_states)-kron(hc1',ones(size(a_grid)))>0)).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc1',ones(size(a_grid))))./omega_bar).^(-gamma)),0)*prob1,...
    max(((((1+r).*repmat(a_grid,1,hc_states)-kron(hc2',ones(size(a_grid))))>0).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc2',ones(size(a_grid))))./omega_bar).^(-gamma)),0)*prob2,...
    max(((((1+r).*repmat(a_grid,1,hc_states)-kron(hc3',ones(size(a_grid))))>0).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc3',ones(size(a_grid))))./omega_bar).^(-gamma)),0)*prob3,...
    max(((((1+r).*repmat(a_grid,1,hc_states)-kron(hc4',ones(size(a_grid))))>0).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc4',ones(size(a_grid))))./omega_bar).^(-gamma)),0)*prob4];


Health_Tran_T= Health_Transition(T-1); %T-1 to T transition
deriv_bequest=[Health_Tran_T(1,:)*deriv_bequest_vec';Health_Tran_T(2,:)*deriv_bequest_vec';Health_Tran_T(3,:)*deriv_bequest_vec']';

% Why is deriv_bequest not smooth?  It has to do with the fact that if you
% carry positive wealth to the next state, in some states government care is are binding,
% and in other states you are not.  Since the derivative in the government
% care state is 0, one would not, as the derivative switches from 0 to
% positive (i.e., the constraint doesn't bind), then there is a jump.  This
% continues until all the jumps are finished, and then it decreases
% (decreasing marginal utility is a functtional assumption.

%% 3c) Apply endogenous grid Method
% Here I generate the endogenous grid.  It comes from inverting the FOC.  I
% do this for each of the 3 states {0,1,2} where the differences emerge
% from the state dependent utility and the differences in continuation
% value.

E_grid_T1=(   kron(((deriv_bequest(:,1)).^(-1/sigma)-y+a_grid),ones(1,hc_states))+kron(hc1_minus',ones(size(a_grid)))   )  /  (1+r);
E_grid_T2=(   kron(((deriv_bequest(:,2)).^(-1/sigma)-y+a_grid),ones(1,hc_states))+kron(hc2_minus',ones(size(a_grid)))   )  /  (1+r);
E_grid_T3=(kron(((deriv_bequest(:,3)./delta).^(-1/sigma)-y+a_grid),ones(1,hc_states))+kron(hc3_minus',ones(size(a_grid)))) /  (1+r);

E_grid_T1(E_grid_T1==Inf)=0*Inf;   E_grid_T2(E_grid_T2==Inf)=0*Inf;   E_grid_T3(E_grid_T3==Inf)=0*Inf;
%If FOC was not invertible because the derivative of the continuation value
%was 0 then I set it to be NaN.  Note that the derivative can be zero
%because of the government care option.

%% Here I find the lower and upper bounds on the non concave regions.
% Any final grid point above concave_i_H  and anything below concave_i_L
% has a one-one mapping from the FOC, so we don't have to check for a
% global max.  Anything in G_NC_i is in a nonconcave part of the
% continuation value, and so there are multiple starting assets that map to
% the same end asset.  We thus need to check for global maxes.
%
for i=1:3
    [concave_L{i}, concave_H{i}]=find_concave(deriv_bequest(:,i));
    G_NC{i}=(concave_L{i}:1:concave_H{i})';
    G_C_L{i}=(1:1:(concave_L{i}-1))';
    G_C_H{i}=((concave_H{i}+1):1:m)';
    store_NC{i}=zeros(length(G_NC{i}),hc_states);
    max_NC{i}=zeros(length(G_NC{i}),hc_states);
    NC_Ind{i}=ismember((1:m)',G_NC{i});
    a_NC{i}=a_grid(G_NC{i}); %Check this line Joseph... Should this be NC_IND? not G_NC? %Checked and it is correct.
%     All_Concave_I{i}=0
    All_Concave_I{i}=0;
    All_Concave_I{i}=(concave_L{i}> concave_H{i});
   
end

% Parallel_Indic=[(1:m*hc_states)',kron((1:hc_states)',ones(m,1))];
NC_Indices_1=cumsum(ismember((1:m)',G_NC{1}));
NC_Indices_2=cumsum(ismember((1:m)',G_NC{2}));
NC_Indices_3=cumsum(ismember((1:m)',G_NC{3}));

keepmat1=zeros(m,3);keepmat2=zeros(m,3);keepmat3=zeros(m,3);
for ii=1:m
    store1=zeros(1,hc_states);store2=zeros(1,hc_states);store3=zeros(1,hc_states);
    max1=zeros(1,hc_states);max2=zeros(1,hc_states);max3=zeros(1,hc_states);
    I1_temp=zeros(1,hc_states); I2_temp=zeros(1,hc_states);I3_temp=zeros(1,hc_states);
    L1=length(G_NC{1});L2=length(G_NC{2}); L3=length(G_NC{3});
    keep1=zeros(1,3);keep2=zeros(1,3);keep3=zeros(1,3);
    for hc_i=1:hc_states
        if (NC_Ind{1}(ii)==1 && All_Concave_I{1}==0)
            
            zi=E_grid_T1(ii,hc_i)*ones(length(G_NC{1}),1);
            c=(1+r)*zi+y-a_NC{1}-hc1_minus(hc_i);
            [~,I]=max(Value_Function1(c, zi, hc1_minus(hc_i)));
            store1(hc_i)=G_NC{1}(I)*(1-isnan(zi(1)));
            max1(hc_i)=a_NC{1}(I);
            Max_Index1=find(ismember(G_NC{1},store1(hc_i)));% This needs to be changed: G_NC{1,jj,kk}
            if (isempty(Max_Index1)==0)
                I1_temp(hc_i)=Max_Index1;
            end
            keep1(hc_i)=(I1_temp(hc_i)==NC_Indices_1(ii)); %Can this be simplified to  keep1(hc_i)=(I==NC_Indices_1(ii));
        end
        
        if (NC_Ind{2}(ii)==1 && All_Concave_I{2}==0)
            
            zi=E_grid_T2(ii,hc_i)*ones(length(G_NC{2}),1);
            c=(1+r)*zi+y-a_NC{2}-hc2_minus(hc_i);
            [~,I]=max(Value_Function2(c, zi, hc2_minus(hc_i)));
            store2(hc_i)=G_NC{2}(I)*(1-isnan(zi(1)));
            max2(hc_i)=a_NC{2}(I);
            Max_Index2=find(ismember(G_NC{2},store2(hc_i))); % Check this... See above
            if (isempty(Max_Index2)==0)
                I2_temp(hc_i)=Max_Index2;
            end
            keep2(hc_i)=(I2_temp(hc_i)==NC_Indices_2(ii));
        end
        
        if (NC_Ind{3}(ii)==1 && All_Concave_I{3}==0)
            zi=E_grid_T3(ii,hc_i)*ones(length(G_NC{3}),1);
            c=(1+r)*zi+y-a_NC{3}-hc3_minus(hc_i);
            [~,I]=max(Value_Function3(c, zi, hc3_minus(hc_i)));
            store3(hc_i)=G_NC{3}(I)*(1-isnan(zi(1)));
            max3(hc_i)=a_NC{3}(I);
            Max_Index3=find(ismember(G_NC{3},store3(hc_i))); %Check this
            if (isempty(Max_Index3)==0)
                I3_temp(hc_i)=Max_Index3;
            end
            keep3(hc_i)=(I3_temp(hc_i)==NC_Indices_3(ii));
        end
        
    end
    keepmat1(ii,:)=keep1;
    keepmat2(ii,:)=keep2;
    keepmat3(ii,:)=keep3;
end
a_grid1=repmat(a_grid,hc_states);    a_grid2=repmat(a_grid,hc_states);    a_grid3=repmat(a_grid,hc_states);

for i=1:(hc_states)
    
    if (All_Concave_I{1}==1)
        I1{i}=(1:m)';
    else
        I1{i}=[G_C_L{1}; find(keepmat1(:,i));G_C_H{1}];
    end
    
    if (All_Concave_I{2}==1)
        I2{i}=(1:m)';
    else
        I2{i}=[G_C_L{2}; find(keepmat2(:,i));G_C_H{2}];
    end
    if (All_Concave_I{3}==1)
        I3{i}=(1:m)';
    else
        I3{i}=[G_C_L{3}; find(keepmat3(:,i));G_C_H{3}];
    end
    
    Endog_1{i}=[E_grid_T1(I1{i},i), a_grid(I1{i})];
    Endog_2{i}=[E_grid_T2(I2{i},i), a_grid(I2{i})];
    Endog_3{i}=[E_grid_T3(I3{i},i), a_grid(I3{i})];
    
    if size(Endog_1{i},1)>1
        %         Policy_store1{i}=interp1_parallel(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear');
        Policy_store1{i}=interp1(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear');
    else
        Policy_store1{i}=a_grid.*NaN;
    end
    
    if size(Endog_2{i},1)>1
        %         Policy_store2{i}=interp1_parallel(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear');
        Policy_store2{i}=interp1(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear');
    else
        Policy_store2{i}=a_grid.*NaN;
    end
    
    
    if size(Endog_3{i},1)>1
        %         Policy_store3-{i}=interp1_parallel(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear');
        Policy_store3{i}=interp1(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear');
    else
        Policy_store3{i}=a_grid.*NaN;
    end
    
    %Here we save which members of our original grid are not included in
    %our endogenous grid.
    I1_skip{i}=ismember(a_grid1(:,i),Endog_1{i}(:,2));
    I2_skip{i}=ismember(a_grid2(:,i),Endog_2{i}(:,2));
    I3_skip{i}=ismember(a_grid3(:,i),Endog_3{i}(:,2));
    
    %Next we store the indices of those grid points that didn't make it to
    %the endogenous grid.
    Miss_1{i}=find(I1_skip{i}==0);
    Miss_2{i}=find(I2_skip{i}==0);
    Miss_3{i}=find(I3_skip{i}==0);
    
    %Here we take the indices of the missing points and find the lower
    %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
    %endogenous grid, but a_grid(Miss_1-1) is, then we save the this point
    %here.  This will be used to define the regions we run a grid search
    %over.
    lower_GS1{i}=ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),Endog_1{i}(:,2)).*a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i);
    lower_GS1{i}(find(-ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),Endog_1{i}(:,2))+1))=[];

                
    %Here we take the indices of the missing points and find the upper
    %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
    %endogenous grid, but a_grid(Miss_1+1) is, then we save the this point
    %here.  This will be used to define the regions we run a grid search
    %over.
    upper_GS1{i}=ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i),Endog_1{i}(:,2)).*a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i);
    upper_GS1{i}(find(-ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1))+1,i),Endog_1{i}(:,2))+1))=[];
    

    %Now we save the the upper and lower bounds of the assets (not indices)
    %that we will run grid searches over.
    lower_GS1A{i}=Endog_1{i}(ismember(Endog_1{i}(:,2),lower_GS1{i}),1);
    upper_GS1A{i}=Endog_1{i}(ismember(Endog_1{i}(:,2),upper_GS1{i}),1);
    
    %The above steps are repeated for states {1,2}
    lower_GS2{i}=ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),Endog_2{i}(:,2)).*a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i);
    lower_GS2{i}(find(-ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),Endog_2{i}(:,2))+1))=[];
    
    upper_GS2{i}=ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),Endog_2{i}(:,2)).*a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i);
    upper_GS2{i}(find(-ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),Endog_2{i}(:,2))+1))=[];
    
    lower_GS2A{i}=Endog_2{i}(ismember(Endog_2{i}(:,2),lower_GS2{i}),1);
    upper_GS2A{i}=Endog_2{i}(ismember(Endog_2{i}(:,2),upper_GS2{i}),1);
    
    lower_GS3{i}=ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),Endog_3{i}(:,2)).*a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i);
    lower_GS3{i}(find(-ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),Endog_3{i}(:,2))+1))=[];
    upper_GS3{i}=ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),Endog_3{i}(:,2)).*a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i);
    upper_GS3{i}(find(-ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),Endog_3{i}(:,2))+1))=[];
    
    lower_GS3A{i}=Endog_3{i}(ismember(Endog_3{i}(:,2),lower_GS3{i}),1);
    upper_GS3A{i}=Endog_3{i}(ismember(Endog_3{i}(:,2),upper_GS3{i}),1);
    %The above steps are repeated for states {1,2}
    % This step adds in the point 0 to the lower bound asset sets if the Endogenous
    % grid we generated before stops at a positive number.  This just means
    % that we will run a grid search over [0,min_endog_grid].
    
    % Finally, we  find all of the assets that lie within the bounds which
    % we will want to run grid searches over.
    
    GS_I1{i}=0;GS_I2{i}=0;GS_I3{i}=0;
    if (isempty(Miss_1{i})==0)
        if (isempty(lower_GS1A{i})==0)
            if (lower_GS1A{i}(end)==Endog_1{i}(end,1) )
                upper_GS1A{i}=[upper_GS1A{i}; a_grid(end)];
            end
        end
        if (isempty(upper_GS1A{i})==0)
            if (upper_GS1A{i}(1)==Endog_1{i}(1,1))
                lower_GS1A{i}=[0; lower_GS1A{i}];
            end
        end
        
        for k=1:length(lower_GS1A{i})
            GS_I1{i}=[GS_I1{i}; a_grid(a_grid>=lower_GS1A{i}(k) & a_grid<=upper_GS1A{i}(k))];
        end
        
    end
    
    if (isempty(Miss_2{i})==0)
        if (isempty(lower_GS2A{i})==0)
            if (lower_GS2A{i}(end)==Endog_2{i}(end,1))
                upper_GS2A{i}=[upper_GS2A{i}; a_grid(end)];
            end
        end
        if (isempty(upper_GS2A{i})==0)
            if (upper_GS2A{i}(1)==Endog_2{i}(1,1))
                lower_GS2A{i}=[0; lower_GS2A{i}];
            end
        end
        
        for k=1:length(lower_GS2A{i})
            GS_I2{i}=[GS_I2{i}; a_grid(a_grid>=lower_GS2A{i}(k) & a_grid<=upper_GS2A{i}(k))];
        end
    end
    
    if (isempty(Miss_3{i})==0)
        if (isempty(lower_GS3A{i})==0)
            if (lower_GS3A{i}(end)==Endog_3{i}(end,1))
                upper_GS3A{i}=[upper_GS3A{i}; a_grid(end)];
            end
        end
        if (isempty(upper_GS3A{i})==0)
            if (upper_GS3A{i}(1)==Endog_3{i}(1,1) )
                lower_GS3A{i}=[0; lower_GS3A{i}];
            end
        end
        
        for k=1:length(lower_GS3A{i})
            GS_I3{i}=[GS_I3{i}; a_grid(a_grid>=lower_GS3A{i}(k) & a_grid<=upper_GS3A{i}(k))];
        end
    end
    
    
    
    %%
    
    GS_I1{i}(1)=[];GS_I2{i}(1)=[];GS_I3{i}(1)=[];
    %%
    
    %Define the budget_bound as the maximum number of assets that one could
    %carry over to the next period.  For asset holdings below the
    %budget_bound we delete from our grid.  We do not need to solve these
    %for the policy functions since the agent has no choice but to go on government care.
    budget_bound{i}=[find(((1+r)*a_grid+y-hc1_minus(i))<0, 1, 'last'), ...
        find(((1+r)*a_grid+y-hc2_minus(i))<0, 1, 'last' ),...
        find(((1+r)*a_grid+y-hc3_minus(i)-chi_LTC)<0, 1, 'last' )];
    GS_I1{i}(GS_I1{i}<a_grid(budget_bound{i}(:,1)))=[];
    GS_I2{i}(GS_I2{i}<a_grid(budget_bound{i}(:,2)))=[];
    GS_I3{i}(GS_I3{i}<a_grid(budget_bound{i}(:,3)))=[];
    
    Out_of_budget1{i}=1:1:(budget_bound{i}(:,1));
    Out_of_budget2{i}=1:1:(budget_bound{i}(:,2));
    Out_of_budget3{i}=1:1:(budget_bound{i}(:,3));
    %     if (isempty(Out_of_budget1{i})==0)
    %         Policy_store1{i}(Out_of_budget1{i})=0; %Set unfeasible points to 0 policy rule
    %     end
    %
    %     if (isempty(Out_of_budget2{i})==0)
    %         Policy_store2{i}(Out_of_budget2{i})=0;
    %     end
    %     if (isempty(Out_of_budget3{i})==0)
    %         Policy_store3{i}(Out_of_budget3{i})=0;
    %     end
    %
end
%%

Policy_joined1=zeros(m,3);Policy_joined2=zeros(m,3);Policy_joined3=zeros(m,3);
for iii=1:m
    asset=a_grid(iii);
    Policy_store_out1=zeros(1,3);Policy_store_out2=zeros(1,3);Policy_store_out3=zeros(1,3);
    for hc_i=1:hc_states
        if  (ismember(asset,GS_I1{hc_i})==1)
            c_gridT=linspace(0,((1+r)*asset+y-hc1_minus(hc_i)),1000)';
            a_start=asset'*ones(length(c_gridT),1);
            [~,point] =max(Value_Function1(c_gridT,a_start,hc1_minus(hc_i)));
            Policy_store_out1(hc_i)=(1+r)*asset+y-hc1_minus(hc_i)-c_gridT(point);
            
        end
        
        
        if  (ismember(asset,GS_I2{hc_i})==1)
            c_gridT=linspace(0,((1+r)*asset+y-hc2_minus(hc_i)),1000)';
            a_start=asset'*ones(length(c_gridT),1);
            [~,point] =max(Value_Function2(c_gridT,a_start,hc2_minus(hc_i)));
            Policy_store_out2(hc_i)=(1+r)*asset+y-hc2_minus(hc_i)-c_gridT(point);
            
        end
        
        if  (ismember(asset,GS_I3{hc_i})==1)
            e_gridT=linspace(0,((1+r)*asset+y-hc3_minus(hc_i)),1000)';
            a_start=asset'*ones(length(e_gridT),1);
            [~,point] =max(Value_Function3(e_gridT,a_start,hc3_minus(hc_i)));
            Policy_store_out3(hc_i)=(1+r)*asset+y-hc3_minus(hc_i)-e_gridT(point);
        end
    end
    Policy_joined1(iii,:)=Policy_store_out1;
    Policy_joined2(iii,:)=Policy_store_out2;
    Policy_joined3(iii,:)=Policy_store_out3;
end
%%




for hc_i=1:hc_states
    Policy_store1{hc_i}(ismember(a_grid,GS_I1{hc_i}))=Policy_joined1(ismember(a_grid,GS_I1{hc_i}),hc_i);
    Policy_store2{hc_i}(ismember(a_grid,GS_I2{hc_i}))=Policy_joined2(ismember(a_grid,GS_I2{hc_i}),hc_i);
    Policy_store3{hc_i}(ismember(a_grid,GS_I3{hc_i}))=Policy_joined3(ismember(a_grid,GS_I3{hc_i}),hc_i);
    
    
    %% Now we calculate the consumption allocations implied by our policy
    % functions.
    Consumption_store1{hc_i}=(1+r)*a_grid+y-hc1_minus(hc_i)-Policy_store1{hc_i};
    Consumption_store2{hc_i}=(1+r)*a_grid+y-hc2_minus(hc_i)-Policy_store2{hc_i};
    Consumption_store3{hc_i}=(1+r)*a_grid+y-hc3_minus(hc_i)-Policy_store3{hc_i};
    if (isempty(Out_of_budget1{hc_i})==0)
        Consumption_store1{hc_i}(Out_of_budget1{hc_i})=0;
    end
    if (isempty(Out_of_budget2{hc_i})==0)
        Consumption_store2{hc_i}(Out_of_budget2{hc_i})=0;
    end
    if (isempty(Out_of_budget3{hc_i})==0)
        Consumption_store3{hc_i}(Out_of_budget3{hc_i})=0;
    end
    
    
    
    % Calculate Value Function values over the grids.
    Value_store1{hc_i}=Value_Function1(Consumption_store1{hc_i}, a_grid, hc1_minus(hc_i));
    Value_store2{hc_i}=Value_Function2(Consumption_store2{hc_i}, a_grid, hc2_minus(hc_i));
    Value_store3{hc_i}=Value_Function3(Consumption_store3{hc_i}, a_grid, hc3_minus(hc_i));
    
    
    
    
    
    %     Check boundary condition: Check and see if the agent is better
    %     off consuming everything or s king with the interior policy
    %     function we calculated previously.
    [Value_store1{hc_i},Max_I1_T{hc_i}]=max([Value_Function1(Consumption_store1{hc_i},a_grid, hc1_minus(hc_i)), ...
        Value_Function1((1+r)*a_grid+y-hc1_minus(hc_i),a_grid, hc1_minus(hc_i))],[],2);
    [Value_store2{hc_i},Max_I2_T{hc_i}]=max([Value_Function2(Consumption_store2{hc_i},a_grid, hc2_minus(hc_i)), ...
        Value_Function2((1+r)*a_grid+y-hc2_minus(hc_i),a_grid, hc2_minus(hc_i))],[],2);
    [Value_store3{hc_i},Max_I3_T{hc_i}]=max([Value_Function3(Consumption_store3{hc_i},a_grid, hc3_minus(hc_i)), ...
        Value_Function3((1+r)*a_grid+y-hc3_minus(hc_i),a_grid, hc3_minus(hc_i))],[],2);
    %Store the points where the agent is better off consuming everything
    Boundary_I1_T=find(Max_I1_T{hc_i}==2); Boundary_I2_T=find(Max_I2_T{hc_i}==2); Boundary_I3_T=find(Max_I3_T{hc_i}==2);
    %and set the policy function of these grid points to zero
    Policy_store1{hc_i}(Boundary_I1_T)=0;
    Policy_store2{hc_i}(Boundary_I2_T)=0;
    Policy_store3{hc_i}(Boundary_I3_T)=0; %and set the policy function of these grid points to zero
    %and the consumption policy of these grid points to all of wealth.
    Consumption_store1{hc_i}(Boundary_I1_T)=(1+r)*a_grid(Boundary_I1_T)+y-hc1_minus(hc_i);
    Consumption_store2{hc_i}(Boundary_I2_T)=(1+r)*a_grid(Boundary_I2_T)+y-hc2_minus(hc_i);
    Consumption_store3{hc_i}(Boundary_I3_T)=(1+r)*a_grid(Boundary_I3_T)+y-hc3_minus(hc_i);
    
    % Check constraint in 3rd health state
    
    Con_fl3{hc_i}=find(Consumption_store3{hc_i}<chi_LTC);
    if (isempty(Out_of_budget3{hc_i})==0)
        Con_fl3{hc_i}(Out_of_budget3{hc_i}<=Con_fl3{hc_i}(end))=[];
    end
    %If the constraint is violated then we run a grid search to see what
    %the optimal policy function is in the constrained consumption set.
    if (isempty(Con_fl3{hc_i})==0)
        for j=Con_fl3{hc_i}'
            asset=a_grid(j);
            e_grid=(chi_LTC:.1:((1+r)*asset+y-hc3_minus(hc_i)))';
            a_start=asset'*ones(length(e_grid),1);
            [~,point] =max(Value_Function3(e_grid,a_start,hc3_minus(hc_i)));
            Policy_store3{hc_i}(j)=(1+r)*asset+y-hc3_minus(hc_i)-e_grid(point);
            Consumption_store3{hc_i}(j)=e_grid(point);
        end
    end
    
    
    
    % We reassign the budget infeasible asset  and consumption policies to
    % 0 and the value function grid to be very small.  These people must go
    % on government care.
    if (isempty(Out_of_budget1{hc_i})==0)
        Policy_store1{hc_i}(Out_of_budget1{hc_i})=0;
        Consumption_store1{hc_i}(Out_of_budget1{hc_i})=0;
        Value_store1{hc_i}(Out_of_budget1{hc_i})=-1e20;
    end
    if (isempty(Out_of_budget2{hc_i})==0)
        Policy_store2{hc_i}(Out_of_budget2{hc_i})=0;
        Consumption_store2{hc_i}(Out_of_budget2{hc_i})=0;
        Value_store2{hc_i}(Out_of_budget2{hc_i})=-1e20;
    end
    if  (isempty(Out_of_budget3{hc_i})==0)
        Policy_store3{hc_i}(Out_of_budget3{hc_i})=0;
        Consumption_store3{hc_i}(Out_of_budget3{hc_i})=0;
        Value_store3{hc_i}(Out_of_budget3{hc_i})=-1e20;
    end
    
    
    
    
    %% 3h) Government Care
    %Here I calculate the value of government care in each state according
    %to the allocations.
    Gov_bequest=[prob1'*bequest(max(-hc1,0),omega_bar,gamma,phi), prob2'*bequest(max(-hc2,0),omega_bar,gamma,phi), prob3'*bequest(max(-hc3,0),omega_bar,gamma,phi) , prob4'*bequest(max(-hc4,0),omega_bar,gamma,phi)]';
    gc1= Utility(C_F, 1, 1,[1 0 0],p)+EYE(1,:)*betta*Health_Transition(T-1)*Gov_bequest;
    gc2= Utility(C_F, 1, 1,[1 0 0],p)+EYE(2,:)*betta*Health_Transition(T-1)*Gov_bequest;
    gc3= Utility(1, LTC_pc, 1,[0 1 0],p)+EYE(3,:)*betta*Health_Transition(T-1)*Gov_bequest;
    
    %Here I Find the points on the Value function grid for which it is
    %optimal for an agent to enter government care.  These points are found
    %by comparing the value function grid to the government care option and
    %selecting the higher.  I also reassign asset and consumption policies,
    % to 0.
    Gov_I1{hc_i}=find(Value_store1{hc_i}<gc1); Value_store1{hc_i}(Gov_I1{hc_i})=gc1;Policy_store1{hc_i}(Gov_I1{hc_i})=0; Consumption_store1{hc_i}(Gov_I1{hc_i})=0;
    Gov_I2{hc_i}=find(Value_store2{hc_i}<gc2); Value_store2{hc_i}(Gov_I2{hc_i})=gc2;Policy_store2{hc_i}(Gov_I2{hc_i})=0; Consumption_store2{hc_i}(Gov_I2{hc_i})=0;
    Gov_I3{hc_i}=find(Value_store3{hc_i}<gc3); Value_store3{hc_i}(Gov_I3{hc_i})=gc3;Policy_store3{hc_i}(Gov_I3{hc_i})=0; Consumption_store3{hc_i}(Gov_I3{hc_i})=0;
    
    %Here I store indexes where constraints are binding.
    Con_bind3{hc_i}=(Consumption_store3{hc_i}==chi_LTC);
    Gov_bind_I1{hc_i}=(Value_store1{hc_i}==gc1);
    Gov_bind_I2{hc_i}=(Value_store2{hc_i}==gc2);
    Gov_bind_I3{hc_i}=(Value_store3{hc_i}==gc3);
    
    %Define Value_function for state 3.  There is no decicision to be made
    %here.
    Value_store4{hc_i}=max(bequest((1+r)*a_grid-hc4_minus(hc_i),omega_bar,gamma,phi), bequest(0,omega_bar,gamma,phi));
    
    
end

for k=1:hc_states
    Value_mat1(:,k)=Value_store1{k};
    Value_mat2(:,k)=Value_store2{k};
    Value_mat3(:,k)=Value_store3{k};
    Value_mat4(:,k)=Value_store4{k};
    Policy_mat1(:,k)=Policy_store1{k};
    Policy_mat2(:,k)=Policy_store2{k};
    Policy_mat3(:,k)=Policy_store3{k};
    Consump_mat1(:,k)=Consumption_store1{k};
    Consump_mat2(:,k)=Consumption_store2{k};
    Consump_mat3(:,k)=Consumption_store3{k};
    
    
    V_sj{1,k,T-Final_age+1}=Value_store1{k};
    V_sj{2,k,T-Final_age+1}=Value_store2{k};
    V_sj{3,k,T-Final_age+1}=Value_store3{k};
    V_sj{4,k,T-Final_age+1}=Value_store4{k};
    
    P_sj{1,k,T-Final_age+1}=Policy_store1{k};
    P_sj{2,k,T-Final_age+1}=Policy_store2{k};
    P_sj{3,k,T-Final_age+1}=Policy_store3{k};
    
    
    C_sj{1,k,T-Final_age+1}=Consumption_store1{k};
    C_sj{2,k,T-Final_age+1}=Consumption_store2{k};
    C_sj{3,k,T-Final_age+1}=Consumption_store3{k};
end






for j=2:(T-Final_age)
     
    for i=1:hc_states
        [hc1(i),prob1(i)]=healthcost(1,T-j+1,i);% Here I generate the health costs and probabilities of each cost for time T-j+1.
        [hc2(i),prob2(i)]=healthcost(2,T-j+1,i);
        [hc3(i),prob3(i)]=healthcost(3,T-j+1,i);
        [hc4(i),prob4(i)]=healthcost(4,T-j+1,i);
        
        
        [hc1_minus(i),prob1_minus(i)]=healthcost(1,T-j,i); % Here I generate the health costs and probabilities of each cost for time T-j
        [hc2_minus(i),prob2_minus(i)]=healthcost(2,T-j,i);
        [hc3_minus(i),prob3_minus(i)]=healthcost(3,T-j,i);
        [hc4_minus(i),prob4_minus(i)]=healthcost(4,T-j,i);
        
        
        
        %Here I calculate the derivatives of each of the continuation  values, depending on if
        %one enters state {0,1,2,3} next period.  The derivations of these expressions are presented
        %in the solution algorithm.  I enforce the relevant constraints
        %that are presented in the write-up, and each is calculated
        %analy ally.
        deriv_1_next(:,i)=(1+r)*Consumption_store1{i}.^(-sigma);  deriv_1_next((Gov_bind_I1{i}'),i)=0;
        deriv_2_next(:,i)=(1+r)*Consumption_store2{i}.^(-sigma);  deriv_2_next((Gov_bind_I2{i}'),i)=0;
        deriv_3_next(:,i)=(1+r)*delta.*Consumption_store3{i}.^(-sigma); deriv_3_next((Gov_bind_I3{i}'),i)=0; deriv_3_next((Con_bind3{i}),i)=0;
        deriv_4_next(:,i)=max((1+r).*((phi+(a_grid-hc4(i))./(omega_bar)).^(-gamma)).*((a_grid-hc4(i))>0),0); deriv_4_next(((a_grid-hc4(i))<0),i)=0;
        
        
    end
    
    Health_Tran= Health_Transition(T-j); %Transition from T-j to T-j+1
    %Here I define the continuation values.  The function takes all
    %arguments and returns the diiscounted expected value of continuing in the next period
    %if you are in state {0,1,2}.  It takes expectation over health
    %transition and health cost.
    Cont_1=@(a)ECV(a, 1, T-j, a_grid,Value_mat1,Value_mat2,Value_mat3,Value_mat4,prob1,prob2,prob3,prob4, Health_Tran, betta);
    Cont_2=@(a)ECV(a, 2, T-j, a_grid,Value_mat1,Value_mat2,Value_mat3,Value_mat4,prob1,prob2,prob3,prob4, Health_Tran, betta);
    Cont_3=@(a)ECV(a, 3, T-j, a_grid,Value_mat1,Value_mat2,Value_mat3,Value_mat4,prob1,prob2,prob3,prob4, Health_Tran, betta);
    %     Defines the value function given todays asset holdings, and health
    %     cost value as a function of asset policy (a_t+1)
    Value_Function1=@(a_next,a,h) Utility((1+r)*a+y-a_next-h, ones(size(a)), ones(size(a)),[1 0 0],p)+Cont_1(a_next)';
    Value_Function2=@(a_next,a,h) Utility((1+r)*a+y-a_next-h, ones(size(a)), ones(size(a)),[1 0 0],p)+Cont_2(a_next)';
    Value_Function3=@(a_next,a,h) Utility(ones(size(a)), (1+r)*a+y-a_next-h, ones(size(a)),[0 1 0],p)+Cont_3(a_next)';
    
    
    
    
    % Take expectation of derivative of continuation values over health
    % cost.  For next periods health state i, deriv_i_next*probi  gives the
    % expected derivative of the continuation value given you are in state
    % i.
    deriv_health_integrated=[deriv_1_next*prob1, deriv_2_next*prob2, deriv_3_next*prob3, deriv_4_next*prob4];
    
    %Takes expectation over health transition.  Gives the
    % expected derivative of the continuation value given you are in state
    % i today.
    deriv_Cont_Val=[Health_Tran(1,:)*deriv_health_integrated';Health_Tran(2,:)*deriv_health_integrated'; Health_Tran(3,:)*deriv_health_integrated']';
    
    
    %% Here I generate the endogenous grid.  It comes from inverting the FOC.  I
    % do this for each of the 3 states {0,1,2} where the differences emerge
    % from the state dependent utility and the differences in continuation
    % value.
    
    E_grid_1=(   kron(((deriv_Cont_Val(:,1)).^(-1/sigma)-y+a_grid),ones(1,hc_states))+kron(hc1_minus',ones(size(a_grid)))   )  /  (1+r);
    E_grid_2=(   kron(((deriv_Cont_Val(:,2)).^(-1/sigma)-y+a_grid),ones(1,hc_states))+kron(hc2_minus',ones(size(a_grid)))   )  /  (1+r);
    E_grid_3=(kron(((deriv_Cont_Val(:,3)./delta).^(-1/sigma)-y+a_grid),ones(1,hc_states))+kron(hc3_minus',ones(size(a_grid)))) /  (1+r);
    
    E_grid_1(E_grid_1==Inf)=0*Inf;   E_grid_2(E_grid_2==Inf)=0*Inf;   E_grid_3(E_grid_3==Inf)=0*Inf;
    
    
    %% Here I find the lower and upper bounds on the non concave regions.
    % Any final grid point above concave_i_H  and anything below concave_i_L
    % has a one-one mapping from the FOC, so we don't have to check for a
    % global max.  Anything in G_NC_i is in a nonconcave part of the
    % continuation value, and so there are multiple starting assets that map to
    % the same end asset.  We thus need to check for global maxes.
    %
    for i=1:3
        [concave_L{i}, concave_H{i}]=find_concave(deriv_Cont_Val(:,i));
        G_NC{i}=(concave_L{i}:1:concave_H{i})';
        G_C_L{i}=(1:1:(concave_L{i}-1))';
        G_C_H{i}=((concave_H{i}+1):1:m)';
        store_NC{i}=zeros(length(G_NC{i}),hc_states); 
        max_NC{i}=zeros(length(G_NC{i}),hc_states);
        NC_Ind{i}=ismember((1:m)',G_NC{i});
        a_NC{i}=a_grid(G_NC{i});
        All_Concave_I{i}=( concave_L{i}> concave_H{i});
        if (isempty(All_Concave_I{i})==1)
            All_Concave_I{i}=0;
        end
        
    end
    
    % Parallel_Indic=[(1:m*hc_states)',kron((1:hc_states)',ones(m,1))];
    NC_Indices_1=cumsum(ismember((1:m)',G_NC{1}));
    NC_Indices_2=cumsum(ismember((1:m)',G_NC{2}));
    NC_Indices_3=cumsum(ismember((1:m)',G_NC{3}));
    
    keepmat1=zeros(m,3);keepmat2=zeros(m,3);keepmat3=zeros(m,3);
   for iv=1:m
        store1=zeros(1,hc_states);store2=zeros(1,hc_states);store3=zeros(1,hc_states);
        max1=zeros(1,hc_states);max2=zeros(1,hc_states);max3=zeros(1,hc_states);
        I1_temp=zeros(1,hc_states); I2_temp=zeros(1,hc_states);I3_temp=zeros(1,hc_states);
        L1=length(G_NC{1}); L2=length(G_NC{2}); L3=length(G_NC{3});
        keep1=zeros(1,3);keep2=zeros(1,3);keep3=zeros(1,3);
        
        for hc_i=1:hc_states
            
            if (NC_Ind{1}(iv)==1 && All_Concave_I{1}==0)
                zi=E_grid_1(iv,hc_i)*ones(length(G_NC{1}),1);
                a=a_grid(G_NC{1});
                [~,I]=max(Value_Function1(a, zi, hc1_minus(hc_i)));
                store1(hc_i)=G_NC{1}(I)*(1-isnan(zi(1)));
                max1(hc_i)=a_NC{1}(I);
                Max_Index1=find(ismember(G_NC{1},store1(hc_i))); % Check G_NC{1} vs. L1
                if (isempty(Max_Index1)==0)
                    I1_temp(hc_i)=Max_Index1;
                end
                keep1(hc_i)=(I1_temp(hc_i)==NC_Indices_1(iv));
            end
            
            if (NC_Ind{2}(iv)==1 && All_Concave_I{2}==0)
                zi=E_grid_2(iv,hc_i)*ones(length(G_NC{2}),1);
                a=a_grid(G_NC{2});
                [~,I]=max(Value_Function2(a, zi, hc2_minus(hc_i)));
                store2(hc_i)=G_NC{2}(I)*(1-isnan(zi(1)));
                max2(hc_i)=a_NC{2}(I);
                Max_Index2=find(ismember(G_NC{2},store2(hc_i)));% Check L2 vs. G_NC{2}
                if (isempty(Max_Index2)==0)
                    I2_temp(hc_i)=Max_Index2;
                end
                keep2(hc_i)=(I2_temp(hc_i)==NC_Indices_2(iv));
                
            end
            if (NC_Ind{3}(iv)==1 && All_Concave_I{3}==0)
                zi=E_grid_3(iv,hc_i)*ones(length(G_NC{3}),1);
                a=a_grid(G_NC{3});
                [~,I]=max(Value_Function3(a, zi, hc2_minus(hc_i)));
                store3(hc_i)=G_NC{3}(I)*(1-isnan(zi(1)));
                max3(hc_i)=a_NC{3}(I);
                Max_Index3=find(ismember(G_NC{3},store3(hc_i)));
                if (isempty(Max_Index3)==0)
                    I3_temp(hc_i)=Max_Index3;
                end
                keep3(hc_i)=(I3_temp(hc_i)==NC_Indices_3(iv));
            end
            
        end
        keepmat1(iv,:)=keep1;
        keepmat2(iv,:)=keep2;
        keepmat3(iv,:)=keep3;
    end
    a_grid1=repmat(a_grid,hc_states);    a_grid2=repmat(a_grid,hc_states);    a_grid3=repmat(a_grid,hc_states);

    for i=1:(hc_states)
        
        if (All_Concave_I{1}==1)
            I1{i}=(1:m)';
        else
            I1{i}=[G_C_L{1}; find(keepmat1(:,i));G_C_H{1}];
        end
        
        if (All_Concave_I{2}==1)
            I2{i}=(1:m)';
        else
            I2{i}=[G_C_L{2}; find(keepmat2(:,i));G_C_H{2}];
        end
        if (All_Concave_I{3}==1)
            I3{i}=(1:m)';
        else
            I3{i}=[G_C_L{3}; find(keepmat3(:,i));G_C_H{3}];
        end
        
        Endog_1{i}=[E_grid_1(I1{i},i), a_grid(I1{i})];
        Endog_2{i}=[E_grid_2(I2{i},i), a_grid(I2{i})];
        Endog_3{i}=[E_grid_3(I3{i},i), a_grid(I3{i})];
        
        if size(Endog_1{i},1)>1
            %         Policy_store1{i}=interp1_parallel(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear');
            Policy_store1{i}=interp1(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear');
        else
            Policy_store1{i}=a_grid.*NaN;
        end
        
        if size(Endog_2{i},1)>1
            %         Policy_store2{i}=interp1_parallel(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear');
            Policy_store2{i}=interp1(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear');
        else
            Policy_store2{i}=a_grid.*NaN;
        end
        
        
        if size(Endog_3{i},1)>1
            %         Policy_store3-{i}=interp1_parallel(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear');
            Policy_store3{i}=interp1(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear');
        else
            Policy_store3{i}=a_grid.*NaN;
        end
        
        %Here we save which members of our original grid are not included in
        %our endogenous grid.
        I1_skip{i}=ismember(a_grid1(:,i),Endog_1{i}(:,2));
        I2_skip{i}=ismember(a_grid2(:,i),Endog_2{i}(:,2));
        I3_skip{i}=ismember(a_grid3(:,i),Endog_3{i}(:,2));
        
        %Next we store the indices of those grid points that didn't make it to
        %the endogenous grid.
        Miss_1{i}=find(I1_skip{i}==0);
        Miss_2{i}=find(I2_skip{i}==0);
        Miss_3{i}=find(I3_skip{i}==0);
        
        %Here we take the indices of the missing points and find the lower
        %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
        %endogenous grid, but a_grid(Miss_1-1) is, then we save the this point
        %here.  This will be used to define the regions we run a grid search
        %over.
        lower_GS1{i}=ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),Endog_1{i}(:,2)).*a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i);
        lower_GS1{i}(find(-ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),Endog_1{i}(:,2))+1))=[];
        
        %Here we take the indices of the missing points and find the upper
        %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
        %endogenous grid, but a_grid(Miss_1+1) is, then we save the this point
        %here.  This will be used to define the regions we run a grid search
        %over.
        upper_GS1{i}=ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i),Endog_1{i}(:,2)).*a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i);
        upper_GS1{i}(find(-ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1))+1,i),Endog_1{i}(:,2))+1))=[];
        
        %Now we save the the upper and lower bounds of the assets (not indices)
        %that we will run grid searches over.
        lower_GS1A{i}=Endog_1{i}(ismember(Endog_1{i}(:,2),lower_GS1{i}),1);
        upper_GS1A{i}=Endog_1{i}(ismember(Endog_1{i}(:,2),upper_GS1{i}),1);
        
        %The above steps are repeated for states {1,2}
        lower_GS2{i}=ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),Endog_2{i}(:,2)).*a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i);
        lower_GS2{i}(find(-ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),Endog_2{i}(:,2))+1))=[];
        
        upper_GS2{i}=ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),Endog_2{i}(:,2)).*a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i);
        upper_GS2{i}(find(-ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),Endog_2{i}(:,2))+1))=[];
        
        lower_GS2A{i}=Endog_2{i}(ismember(Endog_2{i}(:,2),lower_GS2{i}),1);
        upper_GS2A{i}=Endog_2{i}(ismember(Endog_2{i}(:,2),upper_GS2{i}),1);
        
        lower_GS3{i}=ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),Endog_3{i}(:,2)).*a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i);
        lower_GS3{i}(find(-ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),Endog_3{i}(:,2))+1))=[];
        upper_GS3{i}=ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),Endog_3{i}(:,2)).*a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i);
        upper_GS3{i}(find(-ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),Endog_3{i}(:,2))+1))=[];
        
        lower_GS3A{i}=Endog_3{i}(ismember(Endog_3{i}(:,2),lower_GS3{i}),1);
        upper_GS3A{i}=Endog_3{i}(ismember(Endog_3{i}(:,2),upper_GS3{i}),1);
        %The above steps are repeated for states {1,2}
        % This step adds in the point 0 to the lower bound asset sets if the Endogenous
        % grid we generated before stops at a positive number.  This just means
        % that we will run a grid search over [0,min_endog_grid].
        
        % Finally, we  find all of the assets that lie within the bounds which
        % we will want to run grid searches over.
        
        GS_I1{i}=0;GS_I2{i}=0;GS_I3{i}=0;
        if (isempty(Miss_1{i})==0)
            if (isempty(lower_GS1A{i})==0)
                if (lower_GS1A{i}(end)==Endog_1{i}(end,1) )
                    upper_GS1A{i}=[upper_GS1A{i}; a_grid(end)];
                end
            end
            if (isempty(upper_GS1A{i})==0)
                if (upper_GS1A{i}(1)==Endog_1{i}(1,1))
                    lower_GS1A{i}=[0; lower_GS1A{i}];
                end
            end
            
            for k=1:length(lower_GS1A{i})
                GS_I1{i}=[GS_I1{i}; a_grid(a_grid>=lower_GS1A{i}(k) & a_grid<=upper_GS1A{i}(k))];
            end
            
        end
        
        if (isempty(Miss_2{i})==0)
            if (isempty(lower_GS2A{i})==0)
                if (lower_GS2A{i}(end)==Endog_2{i}(end,1))
                    upper_GS2A{i}=[upper_GS2A{i}; a_grid(end)];
                end
            end
            if (isempty(upper_GS2A{i})==0)
                if (upper_GS2A{i}(1)==Endog_2{i}(1,1))
                    lower_GS2A{i}=[0; lower_GS2A{i}];
                end
            end
            
            for k=1:length(lower_GS2A{i})
                GS_I2{i}=[GS_I2{i}; a_grid(a_grid>=lower_GS2A{i}(k) & a_grid<=upper_GS2A{i}(k))];
            end
        end
        
        if (isempty(Miss_3{i})==0)
            if (isempty(lower_GS3A{i})==0)
                if (lower_GS3A{i}(end)==Endog_3{i}(end,1))
                    upper_GS3A{i}=[upper_GS3A{i}; a_grid(end)];
                end
            end
            if (isempty(upper_GS3A{i})==0)
                if (upper_GS3A{i}(1)==Endog_3{i}(1,1) )
                    lower_GS3A{i}=[0; lower_GS3A{i}];
                end
            end
            
            for k=1:length(lower_GS3A{i})
                GS_I3{i}=[GS_I3{i}; a_grid(a_grid>=lower_GS3A{i}(k) & a_grid<=upper_GS3A{i}(k))];
            end
        end
        
        
        
        %%
        
        GS_I1{i}(1)=[];GS_I2{i}(1)=[];GS_I3{i}(1)=[];
        %%
        
        %Define the budget_bound as the maximum number of assets that one could
        %carry over to the next period.  For asset holdings below the
        %budget_bound we delete from our grid.  We do not need to solve these
        %for the policy functions since the agent has no choice but to go on government care.
        budget_bound{i}=[find(((1+r)*a_grid+y-hc1_minus(i))<0,1,'last'), ...
            find(((1+r)*a_grid+y-hc2_minus(i))<0,1,'last'),...
            find(((1+r)*a_grid+y-hc3_minus(i)-chi_LTC)<0,1,'last')];
        GS_I1{i}(GS_I1{i}<a_grid(budget_bound{i}(:,1)))=[];
        GS_I2{i}(GS_I2{i}<a_grid(budget_bound{i}(:,2)))=[];
        GS_I3{i}(GS_I3{i}<a_grid(budget_bound{i}(:,3)))=[];
        
        Out_of_budget1{i}=1:1:(budget_bound{i}(:,1));
        Out_of_budget2{i}=1:1:(budget_bound{i}(:,2));
        Out_of_budget3{i}=1:1:(budget_bound{i}(:,3));
        
    end
    %%
    
    Policy_joined1=zeros(m,3);Policy_joined2=zeros(m,3);Policy_joined3=zeros(m,3);
    for v=1:m
        asset=a_grid(v);
        Policy_store_out1=zeros(1,3);Policy_store_out2=zeros(1,3);Policy_store_out3=zeros(1,3);
        for hc_i=1:hc_states
            if  (ismember(asset,GS_I1{hc_i})==1)
                a_next_grid=linspace(0,((1+r)*asset+y-hc1_minus(hc_i)),1000)';
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function1(a_next_grid,a_start,hc1_minus(hc_i)));
                Policy_store_out1(hc_i)=a_next_grid(point);
                
            end
            
            
            if  (ismember(asset,GS_I2{hc_i})==1)
                a_next_grid=linspace(0,((1+r)*asset+y-hc2_minus(hc_i)),1000)';
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function2(a_next_grid,a_start,hc2_minus(hc_i)));
                Policy_store_out2(hc_i)=a_next_grid(point);
                
            end
            
            if  (ismember(asset,GS_I3{hc_i})==1)
                a_next_grid=linspace(0,((1+r)*asset+y-hc3_minus(hc_i)),1000)';
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function3(a_next_grid,a_start,hc3_minus(hc_i)));
                Policy_store_out3(hc_i)=a_next_grid(point);
            end
        end
        Policy_joined1(v,:)=Policy_store_out1;
        Policy_joined2(v,:)=Policy_store_out2;
        Policy_joined3(v,:)=Policy_store_out3;
    end
    %%
    for hc_i=1:hc_states
        Policy_store1{hc_i}(ismember(a_grid,GS_I1{hc_i}))=Policy_joined1(ismember(a_grid,GS_I1{hc_i}),hc_i);
        Policy_store2{hc_i}(ismember(a_grid,GS_I2{hc_i}))=Policy_joined2(ismember(a_grid,GS_I2{hc_i}),hc_i);
        Policy_store3{hc_i}(ismember(a_grid,GS_I3{hc_i}))=Policy_joined3(ismember(a_grid,GS_I3{hc_i}),hc_i);
        
        
        %% Now we calculate the consumption allocations implied by our policy
        % functions.
        Consumption_store1{hc_i}=(1+r)*a_grid+y-hc1_minus(hc_i)-Policy_store1{hc_i};
        Consumption_store2{hc_i}=(1+r)*a_grid+y-hc2_minus(hc_i)-Policy_store2{hc_i};
        Consumption_store3{hc_i}=(1+r)*a_grid+y-hc3_minus(hc_i)-Policy_store3{hc_i};
        if (isempty(Out_of_budget1{hc_i})==0)
            Consumption_store1{hc_i}(Out_of_budget1{hc_i})=0;
        end
        if (isempty(Out_of_budget2{hc_i})==0)
            Consumption_store2{hc_i}(Out_of_budget2{hc_i})=0;
        end
        if (isempty(Out_of_budget3{hc_i})==0)
            Consumption_store3{hc_i}(Out_of_budget3{hc_i})=0;
        end
        
        
        
        % Calculate Value Function values over the grids.
        Value_store1{hc_i}=Value_Function1(Policy_store1{hc_i}, a_grid, hc1_minus(hc_i));
        Value_store2{hc_i}=Value_Function2(Policy_store2{hc_i}, a_grid, hc2_minus(hc_i));
        Value_store3{hc_i}=Value_Function3(Policy_store3{hc_i}, a_grid, hc3_minus(hc_i));
        
        
        
        
        
        %     Check boundary condition: Check and see if the agent is better
        %     off consuming everything or s king with the interior policy
        %     function we calculated previously.
        [Value_store1{hc_i},Max_I1_T{hc_i}]=max([Value_Function1(Policy_store1{hc_i},a_grid, hc1_minus(hc_i)), ...
            Value_Function1(0,a_grid, hc1_minus(hc_i))],[],2);
        [Value_store2{hc_i},Max_I2_T{hc_i}]=max([Value_Function2(Policy_store2{hc_i},a_grid, hc2_minus(hc_i)), ...
            Value_Function2(0,a_grid, hc2_minus(hc_i))],[],2);
        [Value_store3{hc_i},Max_I3_T{hc_i}]=max([Value_Function3(Policy_store3{hc_i},a_grid, hc3_minus(hc_i)), ...
            Value_Function3(0,a_grid, hc3_minus(hc_i))],[],2);
        %Store the points where the agent is better off consuming everything
        Boundary_I1_T=find(Max_I1_T{hc_i}==2); Boundary_I2_T=find(Max_I2_T{hc_i}==2); Boundary_I3_T=find(Max_I3_T{hc_i}==2);
        %and set the policy function of these grid points to zero
        Policy_store1{hc_i}(Boundary_I1_T)=0;
        Policy_store2{hc_i}(Boundary_I2_T)=0;
        Policy_store3{hc_i}(Boundary_I3_T)=0; %and set the policy function of these grid points to zero
        %and the consumption policy of these grid points to all of wealth.
        Consumption_store1{hc_i}(Boundary_I1_T)=(1+r)*a_grid(Boundary_I1_T)+y-hc1_minus(hc_i);
        Consumption_store2{hc_i}(Boundary_I2_T)=(1+r)*a_grid(Boundary_I2_T)+y-hc2_minus(hc_i);
        Consumption_store3{hc_i}(Boundary_I3_T)=(1+r)*a_grid(Boundary_I3_T)+y-hc3_minus(hc_i);
        
        % Check constraint in 3rd health state
        
        Con_fl3{hc_i}=find(Consumption_store3{hc_i}<chi_LTC);
        if (isempty(Out_of_budget3{hc_i})==0)
            Con_fl3{hc_i}(Out_of_budget3{hc_i}<=max(Con_fl3{hc_i}))=[];
        end
        %If the constraint is violated then we run a grid search to see what
        %the optimal policy function is in the constrained consumption set.
        if (isempty(Con_fl3{hc_i})==0)
            for k=Con_fl3{hc_i}'
                asset=a_grid(k);
                a_next_grid=[(0:.1:((1+r)*asset+y-hc3_minus(i)-chi_LTC))';((1+r)*asset+y-hc3_minus(i)-chi_LTC)];
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function3(a_next_grid,a_start,hc3_minus(hc_i)));
                Policy_store3{hc_i}(k)=a_next_grid(point);
                Consumption_store3{hc_i}(k)=(1+r)*asset+y-hc3_minus(hc_i)-a_next_grid(point);
            end
        end
        
        
        
        % We reassign the budget infeasible asset  and consumption policies to
        % 0 and the value function grid to be very small.  These people must go
        % on government care.
        if (isempty(Out_of_budget1{hc_i})==0)
            Policy_store1{hc_i}(Out_of_budget1{hc_i})=0;
            Consumption_store1{hc_i}(Out_of_budget1{hc_i})=0;
            Value_store1{hc_i}(Out_of_budget1{hc_i})=-1e20;
        end
        if (isempty(Out_of_budget2{hc_i})==0)
            Policy_store2{hc_i}(Out_of_budget2{hc_i})=0;
            Consumption_store2{hc_i}(Out_of_budget2{hc_i})=0;
            Value_store2{hc_i}(Out_of_budget2{hc_i})=-1e20;
        end
        if  (isempty(Out_of_budget3{hc_i})==0)
            Policy_store3{hc_i}(Out_of_budget3{hc_i})=0;
            Consumption_store3{hc_i}(Out_of_budget3{hc_i})=0;
            Value_store3{hc_i}(Out_of_budget3{hc_i})=-1e20;
        end
        
        
        
        
        %% 3h) Government Care
        %Here I calculate the value of government care in each state according
        %to the allocations.
        
        gc1= Utility(C_F, 1, 1,[1 0 0],p)+Cont_1(0);
        gc2= Utility(C_F, 1, 1,[1 0 0],p)+Cont_2(0);
        gc3= Utility(1, LTC_pc, 1,[0 1 0],p)+Cont_3(0);
        
        %Here I Find the points on the Value function grid for which it is
        %optimal for an agent to enter government care.  These points are found
        %by comparing the value function grid to the government care option and
        %selecting the higher.  I also reassign asset and consumption policies,
        % to 0.
        Gov_I1{hc_i}=find(Value_store1{hc_i}<gc1); Value_store1{hc_i}(Gov_I1{i})=gc1;Policy_store1{hc_i}(Gov_I1{i})=0; Consumption_store1{i}(Gov_I1{hc_i})=0;
        Gov_I2{hc_i}=find(Value_store2{hc_i}<gc2); Value_store2{hc_i}(Gov_I2{i})=gc2;Policy_store2{hc_i}(Gov_I2{i})=0; Consumption_store2{i}(Gov_I2{hc_i})=0;
        Gov_I3{hc_i}=find(Value_store3{hc_i}<gc3); Value_store3{hc_i}(Gov_I3{i})=gc3;Policy_store3{hc_i}(Gov_I3{i})=0; Consumption_store3{i}(Gov_I3{hc_i})=0;
        
        %Here I store indexes where constraints are binding.
        Con_bind3{hc_i}=(Consumption_store3{hc_i}==chi_LTC);
        Gov_bind_I1{hc_i}=(Value_store1{hc_i}==gc1);
        Gov_bind_I2{hc_i}=(Value_store2{hc_i}==gc2);
        Gov_bind_I3{hc_i}=(Value_store3{hc_i}==gc3);
        
        %Define Value_function for state 3.  There is no decicision to be made
        %here.
        Value_store4{hc_i}=max(bequest((1+r)*a_grid-hc4_minus(hc_i),omega_bar,gamma,phi), bequest(0,omega_bar,gamma,phi));
        
        
    end
    
    for k=1:hc_states
        Value_mat1(:,k)=Value_store1{k};
        Value_mat2(:,k)=Value_store2{k};
        Value_mat3(:,k)=Value_store3{k};
        Value_mat4(:,k)=Value_store4{k};
        Policy_mat1(:,k)=Policy_store1{k};
        Policy_mat2(:,k)=Policy_store2{k};
        Policy_mat3(:,k)=Policy_store3{k};
        Consump_mat1(:,k)=Consumption_store1{k};
        Consump_mat2(:,k)=Consumption_store2{k};
        Consump_mat3(:,k)=Consumption_store3{k};
        
        
        V_sj{1,k,T-Final_age+2-j}=Value_store1{k};
        V_sj{2,k,T-Final_age+2-j}=Value_store2{k};
        V_sj{3,k,T-Final_age+2-j}=Value_store3{k};
        V_sj{4,k,T-Final_age+2-j}=Value_store4{k};
        
        P_sj{1,k,T-Final_age+2-j}=Policy_store1{k};
        P_sj{2,k,T-Final_age+2-j}=Policy_store2{k};
        P_sj{3,k,T-Final_age+2-j}=Policy_store3{k};
        
        
        C_sj{1,k,T-Final_age+2-j}=Consumption_store1{k};
        C_sj{2,k,T-Final_age+2-j}=Consumption_store2{k};
        C_sj{3,i,T-Final_age+2-j}=Consumption_store3{i};
    end
     
    
   
    breakpoint=1;
end