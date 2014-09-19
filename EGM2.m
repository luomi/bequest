function [V_sj, P_sj,C_sj]=EGM(j,y,T,parameter,r,sigma,betta_in,delta,omega_bar,gamma,phi)
% This function takes asset values a, health state s, age j, and income y,
% and returns the Value function and Policy function (asset holdings for
% period j+1) in these states.

%% 1 - Initialization of parameters
% clc
disp(parameter)    
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

% omega_bar=parameter(4);
% gamma= parameter(5);
% phi= parameter(6);
% 

Final_age=j;
% final_a=a;
% final_h=h;
% final_s=s;


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
% deriv_bequest_vec=[(((((1+r).*repmat(a_grid,1,hc_states)-kron(hc1',ones(size(a_grid)))>0)).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc1',ones(size(a_grid))))./omega_bar).^(-gamma)))*prob1,...
%     (((((1+r).*repmat(a_grid,1,hc_states)-kron(hc2',ones(size(a_grid))))>0).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc2',ones(size(a_grid))))./omega_bar).^(-gamma)))*prob2,...
%     (((((1+r).*repmat(a_grid,1,hc_states)-kron(hc3',ones(size(a_grid))))>0).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc3',ones(size(a_grid))))./omega_bar).^(-gamma)))*prob3,...
%     (((((1+r).*repmat(a_grid,1,hc_states)-kron(hc4',ones(size(a_grid))))>0).*betta.*(1+r).*(phi+((1+r).*repmat(a_grid,1,hc_states)-kron(hc4',ones(size(a_grid))))./omega_bar).^(-gamma)))*prob4];
% 



Health_Tran_T= Health_Transition(T-1);
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
[concave_1_L, concave_1_H]=find_concave(deriv_bequest(:,1));
[concave_2_L, concave_2_H]=find_concave(deriv_bequest(:,2));
[concave_3_L, concave_3_H]=find_concave(deriv_bequest(:,3));

%
G_NC_1=(concave_1_L:1:concave_1_H)';
G_NC_2=(concave_2_L:1:concave_2_H)';
G_NC_3=(concave_3_L:1:concave_3_H)';
G_C_1L=(1:1:(concave_1_L-1))';
G_C_1H=((concave_1_H+1):1:m)';
G_C_2L=(1:1:(concave_2_L-1))';
G_C_2H=((concave_2_H+1):1:m)';
G_C_3L=(1:1:(concave_3_L-1))';
G_C_3H=((concave_3_H+1):1:m)';



%
store1=zeros(length(G_NC_1),hc_states); max1=zeros(length(G_NC_1),hc_states);
store2=zeros(length(G_NC_2),hc_states); max2=zeros(length(G_NC_2),hc_states);
store3=zeros(length(G_NC_3),hc_states); max3=zeros(length(G_NC_3),hc_states);

%% Here we check for each point in the nonconcave section of the grid
% whether for each pair (a_start, a_end) if a_end is the global max out of
% the set (a_end1, a_end2,....) for each a_start.   If the pair does
% constitute a global max pair, then we save it.  (we do this when we check
% whether storei=G_NC_i).  If the pair doesn't constitute a global max pair
% then we discard it.
count=0;
for i=G_NC_1'
    count=count+1;
    a=a_grid(G_NC_1);
    for j=1:hc_states
        zi=E_grid_T1(i,j)*ones(length(G_NC_1),1);
        c=(1+r)*zi+y-a-hc1_minus(j);
        [~,I]=max(Value_Function1(c, zi, hc1_minus(j)));
        store1(count,j)=G_NC_1(I)*(1-isnan(zi(1)));
        max1(count,j)=a(I);
    end
end
[I1_i,I1_j]=find(store1==repmat(G_NC_1,1,hc_states));


count=0;
for i=G_NC_2'
    count=count+1;
    a=a_grid(G_NC_2);
    for j=1:hc_states
        zi=E_grid_T2(i,j)*ones(length(G_NC_2),1);
        c=(1+r)*zi+y-a-hc2_minus(j);
        [~,I]=max(Value_Function2(c, zi, hc2_minus(j)));
        store2(count,j)=G_NC_2(I)*(1-isnan(zi(1)));
        max2(count,j)=a(I);
    end
    
end
[I2_i,I2_j]=find(store2==repmat(G_NC_2,1,hc_states));


count=0;
for i=G_NC_3'
    count=count+1;
    a=a_grid(G_NC_3);
    for j=1:hc_states
        zi=E_grid_T3(i,j)*ones(length(G_NC_3),1);
        c=(1+r)*zi+y-a-hc3_minus(j);
        [~,I]=max(Value_Function3(c, zi, hc3_minus(j)));
        store3(count,j)=G_NC_3(I)*(1-isnan(zi(1)));
        max3(count,j)=a(I);
    end
end
[I3_i,I3_j]=find(store3==repmat(G_NC_3,1,hc_states));
 


%% Here we save the set of points Ii that we've decided to keep.  These
% consist of the concave regions (G_C_iH and G_C_iL) and the points that
% constituted the global max pairs.  Finally, if the whole problem is
% concave,.i.e.  concave_i_L>concave_i_H then we store the whole grid as
% the points we keep.  These are the points we will interpolate over.
for i=1:hc_states
    I1{i}=[G_C_1L; I1_i(I1_j==i); G_C_1H];
    I2{i}=[G_C_2L; I2_i(I2_j==i); G_C_2H];
    I3{i}=[G_C_3L; I3_i(I3_j==i); G_C_3H];
end
if concave_1_L>concave_1_H
    for i=1:hc_states
        I1{i}=(1:m);
    end
end

if concave_2_L>concave_2_H
    for i=1:hc_states
        I2{i}=(1:m);
    end
end

if concave_1_L>concave_1_H
    for i=1:hc_states
        I3{i}=(1:m);
    end
end

%% 3d) Check and make sure that the upper limit of the grid is covered
% If income (y) is large, than an agent may actually be a saver and end up
% with more final period assets than he has at the beginning.  This isn't a
% problem in and of itself, but we need to be careful about grids in this
% case to ensure the upper part of the grid is covered for interpolation
% purposes.  Here, I check if the maximum of our asset grid (max(a_grid))
% is less than the endogenous grid we generated.  If not, then we run a
% grid search to add in this point and ensure coverage.  I do this for each
% of the health states.

a_grid1=repmat(a_grid,1,hc_states); a_grid2=repmat(a_grid,1,hc_states); a_grid3=repmat(a_grid,1,hc_states);


I1_max=(max(E_grid_T1)<max(a_grid)); I2_max=(max(E_grid_T2)<max(a_grid)); I3_max=(max(E_grid_T3)<max(a_grid));

if sum(I1_max)>0
    
    final_asset1=zeros(1,hc_states);
    maxbudget=(1+r)*max(a_grid)+y-hc1_minus;
    
    for i=1:hc_states
        c_grid=0:10:maxbudget(i);
        asset=max(a_grid)*ones(length(c_grid),1);
        [~,point] =max(Value_Function1(c_grid',asset,hc1_minus(i)));
        final_asset1(i)=(1+r)*asset(1)+y-hc1_minus(i)-c_grid(point);
        I1{i}=[I1{i};m+1];
    end
    
    E_grid_T1=[E_grid_T1; max(a_grid)*ones(1,hc_states)];
    a_grid1=[a_grid1; final_asset1];
end


if sum(I2_max)>0
    
    final_asset2=zeros(1,hc_states);
    maxbudget=(1+r)*max(a_grid)+y-hc2_minus;
    
    for i=1:hc_states
        c_grid=0:10:maxbudget(i);
        asset=max(a_grid)*ones(length(c_grid),1);
        [~,point] =max(Value_Function2(c_grid',asset,hc2_minus(i)));
        final_asset2(i)=(1+r)*asset(1)+y-hc2_minus(i)-c_grid(point);
        I2{i}=[I2{i};m+1];
    end
    E_grid_T2=[E_grid_T2; max(a_grid)*ones(1,hc_states)];
    a_grid2=[a_grid2; final_asset2];
end

if sum(I3_max)>0
    
    final_asset3=zeros(1,hc_states);
    maxbudget=(1+r)*max(a_grid)+y-hc3_minus;
    
    for i=1:hc_states
        e_grid=0:10:maxbudget(i);
        asset=max(a_grid)*ones(length(e_grid),1);
        [~,point] =max(Value_Function3(e_grid',asset,hc3_minus(i)));
        final_asset3(i)=(1+r)*asset(1)+y-hc3_minus(i)-e_grid(point);
        I3{i}=[I3{i};m+1];
    end
    
    E_grid_T3=[E_grid_T3; max(a_grid)*ones(1,hc_states)];
    a_grid3=[a_grid3; final_asset3];
end

%% Interpolation
%Because as outlined in the solution algorithm there exist discrete jumps
%in the policy functions.  To prevent interpolating over these jumps, we
%thus only want to interpolate over those pairs that are self-insuring for
%the same things.  As outlined in the solution algorithm, these sets
%correspond to the connected regions of our endogenous grid, i.e. the
%regions where we didn't throw out any points between points.  This step is
%implemented here.

%Here I save the endogenous grid pairs we kept from the previous steps.
for i=1:hc_states
    Endog_1{i}=[E_grid_T1(I1{i},i), a_grid1(I1{i},i)];
    Endog_2{i}=[E_grid_T2(I2{i},i), a_grid2(I2{i},i)];
    Endog_3{i}=[E_grid_T3(I3{i},i), a_grid3(I3{i},i)];
end

%I initially interpolate over all of these pairs.
for i=1:hc_states
    if size(Endog_1{i},1)>1
        Policy_store1{i}=interp1(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear');
    else
        Policy_store1{i}=a_grid.*NaN;
    end
    
    if size(Endog_2{i},1)>1    
        Policy_store2{i}=interp1(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear');
    else
        Policy_store2{i}=a_grid.*NaN;
    end
    
        
    if size(Endog_3{i},1)>1     
        Policy_store3{i}=interp1(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear');
    else
        Policy_store3{i}=a_grid.*NaN;
    end
    
end


for i=1:hc_states
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
    
    
    %Define the budget_bound as the maximum number of assets that one could
    %carry over to the next period.  For asset holdings below the
    %budget_bound we delete from our grid.  We do not need to solve these
    %for the policy functions since the agent has no choice but to go on government care.
    budget_bound{i}=[find(((1+r)*a_grid+y-hc1_minus(i))>0,1), ...
        find(((1+r)*a_grid+y-hc2_minus(i))>0, 1 ),...
        find(((1+r)*a_grid+y-hc3_minus(i)-chi_LTC)>0, 1 )];
    GS_I1{i}(GS_I1{i}<a_grid(budget_bound{i}(:,1)))=[];
    GS_I2{i}(GS_I2{i}<a_grid(budget_bound{i}(:,2)))=[];
    GS_I3{i}(GS_I3{i}<a_grid(budget_bound{i}(:,3)))=[];
    
    Out_of_budget1{i}=1:1:(budget_bound{i}(:,1)-1);
    Out_of_budget2{i}=1:1:(budget_bound{i}(:,2)-1);
    Out_of_budget3{i}=1:1:(budget_bound{i}(:,3)-1);
    if (isempty(Out_of_budget1{i})==0)
        Policy_store1{i}(Out_of_budget1{i})=0; %Set unfeasible points to 0 policy rule
    end
    
    if (isempty(Out_of_budget2{i})==0)
        Policy_store2{i}(Out_of_budget2{i})=0;
    end
    if (isempty(Out_of_budget3{i})==0)
        Policy_store3{i}(Out_of_budget3{i})=0;
        
    end
    % Here we conduct a grid search in the regions that aren't connected
    % and store the polciy functions.
    if (isempty(GS_I1{i})==0)
        for asset=GS_I1{i}'
            c_grid=(0:.1:((1+r)*asset+y-hc1_minus(i)))';
            a_start=asset'*ones(length(c_grid),1);
            [~,point] =max(Value_Function1(c_grid,a_start,hc1_minus(i)));
            Policy_store1{i}(a_grid==asset)=(1+r)*asset+y-hc1_minus(i)-c_grid(point);
            
        end
    end
    
    if (isempty(GS_I2{i})==0)
        for asset=GS_I2{i}'
            c_grid=(0:.1:((1+r)*asset+y-hc2_minus(i)))';
            a_start=asset'*ones(length(c_grid),1);
            [~,point] =max(Value_Function2(c_grid,a_start,hc2_minus(i)));
            Policy_store2{i}(a_grid==asset)=(1+r)*asset+y-hc2_minus(i)-c_grid(point);
            
        end
    end
    
    if (isempty(GS_I3{i})==0)
        for asset=GS_I3{i}'
            e_grid=(0:.1:((1+r)*asset+y-hc3_minus(i)))';
            a_start=asset'*ones(length(e_grid),1);
            [~,point] =max(Value_Function3(e_grid,a_start,hc3_minus(i)));
            Policy_store3{i}(a_grid==asset)=(1+r)*asset+y-hc3_minus(i)-e_grid(point);
        end
    end
    
    
    %% Now we calculate the consumption allocations implied by our policy
    % functions.
    Consumption_store1{i}=(1+r)*a_grid+y-hc1_minus(i)-Policy_store1{i};
    Consumption_store2{i}=(1+r)*a_grid+y-hc2_minus(i)-Policy_store2{i};
    Consumption_store3{i}=(1+r)*a_grid+y-hc3_minus(i)-Policy_store3{i};
    if (isempty(Out_of_budget1{i})==0)
        Consumption_store1{i}(Out_of_budget1{i})=0;
    end
    if (isempty(Out_of_budget2{i})==0)
        Consumption_store2{i}(Out_of_budget2{i})=0;
    end
    if (isempty(Out_of_budget3{i})==0)
        Consumption_store3{i}(Out_of_budget3{i})=0;
    end
    
    
    
    % Calculate Value Function values over the grids.
    Value_store1{i}=Value_Function1(Consumption_store1{i}, a_grid, hc1_minus(i));
    Value_store2{i}=Value_Function2(Consumption_store2{i}, a_grid, hc2_minus(i));
    Value_store3{i}=Value_Function3(Consumption_store3{i}, a_grid, hc3_minus(i));
    
    
    
    
    
    %     Check boundary condition: Check and see if the agent is better
    %     off consuming everything or sticking with the interior policy
    %     function we calculated previously.
    [Value_store1{i},Max_I1_T{i}]=max([Value_Function1(Consumption_store1{i},a_grid, hc1_minus(i)), ...
        Value_Function1((1+r)*a_grid+y-hc1_minus(i),a_grid, hc1_minus(i))],[],2);
    [Value_store2{i},Max_I2_T{i}]=max([Value_Function2(Consumption_store2{i},a_grid, hc2_minus(i)), ...
        Value_Function2((1+r)*a_grid+y-hc2_minus(i),a_grid, hc2_minus(i))],[],2);
    [Value_store3{i},Max_I3_T{i}]=max([Value_Function3(Consumption_store3{i},a_grid, hc3_minus(i)), ...
        Value_Function3((1+r)*a_grid+y-hc3_minus(i),a_grid, hc3_minus(i))],[],2);
    %Store the points where the agent is better off consuming everything
    Boundary_I1_T=find(Max_I1_T{i}==2); Boundary_I2_T=find(Max_I2_T{i}==2); Boundary_I3_T=find(Max_I3_T{i}==2);
    %and set the policy function of these grid points to zero
    Policy_store1{i}(Boundary_I1_T)=0;
    Policy_store2{i}(Boundary_I2_T)=0;
    Policy_store3{i}(Boundary_I3_T)=0; %and set the policy function of these grid points to zero
    %and the consumption policy of these grid points to all of wealth.
    Consumption_store1{i}(Boundary_I1_T)=(1+r)*a_grid(Boundary_I1_T)+y-hc1_minus(i);
    Consumption_store2{i}(Boundary_I2_T)=(1+r)*a_grid(Boundary_I2_T)+y-hc2_minus(i);
    Consumption_store3{i}(Boundary_I3_T)=(1+r)*a_grid(Boundary_I3_T)+y-hc3_minus(i);
    
    % Check constraint in 3rd health state
    
    Con_fl3{i}=find(Consumption_store3{i}<chi_LTC);
    if (isempty(Out_of_budget3{i})==0)
        Con_fl3{i}(Out_of_budget3{i}<=Con_fl3{i}(end))=[];
    end
    %If the constraint is violated then we run a grid search to see what
    %the optimal policy function is in the constrained consumption set.
    if (isempty(Con_fl3{i})==0)
        for j=Con_fl3{i}'
            asset=a_grid(j);
            e_grid=(chi_LTC:.1:((1+r)*asset+y-hc3_minus(i)))';
            a_start=asset'*ones(length(e_grid),1);
            [~,point] =max(Value_Function3(e_grid,a_start,hc3_minus(i)));
            Policy_store3{i}(j)=(1+r)*asset+y-hc3_minus(i)-e_grid(point);
            Consumption_store3{i}(j)=e_grid(point);
        end
    end
    
    
    
    % We reassign the budget infeasible asset  and consumption policies to
    % 0 and the value function grid to be very small.  These people must go
    % on government care.
    if (isempty(Out_of_budget1{i})==0)
        Policy_store1{i}(Out_of_budget1{i})=0;
        Consumption_store1{i}(Out_of_budget1{i})=0;
        Value_store1{i}(Out_of_budget1{i})=-1e20;
    end
    if (isempty(Out_of_budget2{i})==0)
        Policy_store2{i}(Out_of_budget2{i})=0;
        Consumption_store2{i}(Out_of_budget2{i})=0;
        Value_store2{i}(Out_of_budget2{i})=-1e20;
    end
    if  (isempty(Out_of_budget3{i})==0)
        Policy_store3{i}(Out_of_budget3{i})=0;
        Consumption_store3{i}(Out_of_budget3{i})=0;
        Value_store3{i}(Out_of_budget3{i})=-1e20;
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
    Gov_I1{i}=find(Value_store1{i}<gc1); Value_store1{i}(Gov_I1{i})=gc1;Policy_store1{i}(Gov_I1{i})=0; Consumption_store1{i}(Gov_I1{i})=0;
    Gov_I2{i}=find(Value_store2{i}<gc2); Value_store2{i}(Gov_I2{i})=gc2;Policy_store2{i}(Gov_I2{i})=0; Consumption_store2{i}(Gov_I2{i})=0;
    Gov_I3{i}=find(Value_store3{i}<gc3); Value_store3{i}(Gov_I3{i})=gc3;Policy_store3{i}(Gov_I3{i})=0; Consumption_store3{i}(Gov_I3{i})=0;
    
    %Here I store indexes where constraints are binding.
    Con_bind3{i}=(Consumption_store3{i}==chi_LTC);
    Gov_bind_I1{i}=(Value_store1{i}==gc1);
    Gov_bind_I2{i}=(Value_store2{i}==gc2);
    Gov_bind_I3{i}=(Value_store3{i}==gc3);
    
    %Define Value_function for state 3.  There is no decicision to be made
    %here.
    Value_store4{i}=max(bequest((1+r)*a_grid-hc4_minus(i),omega_bar,gamma,phi), bequest(0,omega_bar,gamma,phi));
    
    
end



%This changes the cell-arrays to matrices.  This is done so that we can
%multiply in the next steps to make it easier to calculate expectations.
for i=1:hc_states
    Value_mat1(:,i)=Value_store1{i};
    Value_mat2(:,i)=Value_store2{i};
    Value_mat3(:,i)=Value_store3{i};
    Value_mat4(:,i)=Value_store4{i};
    Policy_mat1(:,i)=Policy_store1{i};
    Policy_mat2(:,i)=Policy_store2{i};
    Policy_mat3(:,i)=Policy_store3{i};
    Consump_mat1(:,i)=Consumption_store1{i};
    Consump_mat2(:,i)=Consumption_store2{i};
    Consump_mat3(:,i)=Consumption_store3{i};
    
    
    V_sj{1,i,T-Final_age+1}=Value_store1{i};
    V_sj{2,i,T-Final_age+1}=Value_store2{i};
    V_sj{3,i,T-Final_age+1}=Value_store3{i};
    V_sj{4,i,T-Final_age+1}=Value_store4{i};
    
    P_sj{1,i,T-Final_age+1}=Policy_store1{i};
    P_sj{2,i,T-Final_age+1}=Policy_store2{i};
    P_sj{3,i,T-Final_age+1}=Policy_store3{i};
    
    
    C_sj{1,i,T-Final_age+1}=Consumption_store1{i};
    C_sj{2,i,T-Final_age+1}=Consumption_store2{i};
    C_sj{3,i,T-Final_age+1}=Consumption_store3{i};
end 


    
%% Now, we solve the T-2 to t problem.
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
        %analytically.
        deriv_1_next(:,i)=(1+r)*Consumption_store1{i}.^(-sigma);  deriv_1_next((Gov_bind_I1{i}'),i)=0;
        deriv_2_next(:,i)=(1+r)*Consumption_store2{i}.^(-sigma);  deriv_2_next((Gov_bind_I2{i}'),i)=0;
        deriv_3_next(:,i)=(1+r)*delta.*Consumption_store3{i}.^(-sigma); deriv_3_next((Gov_bind_I3{i}'),i)=0; deriv_3_next((Con_bind3{i}),i)=0;
        deriv_4_next(:,i)=max((1+r).*((phi+(a_grid-hc4(i))./(omega_bar)).^(-gamma)).*((a_grid-hc4(i))>0),0); deriv_4_next(((a_grid-hc4(i))<0),i)=0;
        
        
    end
    
    Health_Tran= Health_Transition(T-j);
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
    % the same end asset.  We thus need to check for global maxes for these
    % points.
    [concave_1_L, concave_1_H]=find_concave(deriv_Cont_Val(:,1));
    [concave_2_L, concave_2_H]=find_concave(deriv_Cont_Val(:,2));
    [concave_3_L, concave_3_H]=find_concave(deriv_Cont_Val(:,3));
    
    
    
    %
    G_NC_1=(concave_1_L:1:concave_1_H)';
    G_NC_2=(concave_2_L:1:concave_2_H)';
    G_NC_3=(concave_3_L:1:concave_3_H)';
    G_C_1L=(1:1:(concave_1_L-1))';
    G_C_1H=((concave_1_H+1):1:m)';
    G_C_2L=(1:1:(concave_2_L-1))';
    G_C_2H=((concave_2_H+1):1:m)';
    G_C_3L=(1:1:(concave_3_L-1))';
    G_C_3H=((concave_3_H+1):1:m)';
    
    
    
    %
    store1=zeros(length(G_NC_1),hc_states); max1=zeros(length(G_NC_1),hc_states);
    store2=zeros(length(G_NC_2),hc_states); max2=zeros(length(G_NC_2),hc_states);
    store3=zeros(length(G_NC_3),hc_states); max3=zeros(length(G_NC_3),hc_states);
    %% Here we check for each point in the nonconcave section of the grid
    % whether for each pair (a_start, a_end) if a_end is the global max out of
    % the set (a_end1, a_end2,....) for each a_start.   If the pair does
    % constitute a global max pair, then we save it.  (we do this when we check
    % whether storei=G_NC_i).  If the pair doesn't constitute a global max pair
    % then we discard it.
   
    count=0;
    for i=G_NC_1'
        count=count+1;
        a=a_grid(G_NC_1);
        for k=1:hc_states
            zi=E_grid_1(i,k)*ones(length(G_NC_1),1);
            [~,I]=max(Value_Function1(a, zi, hc1_minus(k)));
            store1(count,k)=G_NC_1(I)*(1-isnan(zi(1)));
            max1(count,k)=a(I);
        end
    end
    [I1_i,I1_j]=find(store1==repmat(G_NC_1,1,hc_states));
    
    
    count=0;
    for i=G_NC_2'
        count=count+1;
        a=a_grid(G_NC_2);
        for k=1:hc_states
            zi=E_grid_2(i,k)*ones(length(G_NC_2),1);
            [~,I]=max(Value_Function2(a, zi, hc2_minus(k)));
            store2(count,k)=G_NC_2(I)*(1-isnan(zi(1)));
            max2(count,k)=a(I);
        end
        
    end
    [I2_i,I2_j]=find(store2==repmat(G_NC_2,1,hc_states));
    
    
    count=0;
    for i=G_NC_3'
        count=count+1;
        a=a_grid(G_NC_3);
        for k=1:hc_states
            zi=E_grid_3(i,k)*ones(length(G_NC_3),1);
            [~,I]=max(Value_Function3(a, zi, hc3_minus(k)));
            store3(count,k)=G_NC_3(I)*(1-isnan(zi(1)));
            max3(count,k)=a(I);
        end
    end
    [I3_i,I3_j]=find(store3==repmat(G_NC_3,1,hc_states));
    
    %% Here we save the set of points Ii that we've decided to keep.  These
    % consist of the concave regions (G_C_iH and G_C_iL) and the points that
    % constituted the global max pairs.  Finally, if the whole problem is
    % concave,.i.e.  concave_i_L>concave_i_H then we store the whole grid as
    % the points we keep.  These are the points we will interpolate over.
    
    for i=1:hc_states
        I1{i}=[G_C_1L; G_NC_1(I1_i(I1_j==i)); G_C_1H];
        I2{i}=[G_C_2L; G_NC_2(I2_i(I2_j==i)); G_C_2H];
        I3{i}=[G_C_3L; G_NC_3(I3_i(I3_j==i)); G_C_3H];
    end
    if concave_1_L>concave_1_H
        for i=1:hc_states
            I1{i}=(1:m);
        end
    end
    
    if concave_2_L>concave_2_H
        for i=1:hc_states
            I2{i}=(1:m);
        end
    end
    
    if concave_3_L>concave_3_H
        for i=1:hc_states
            I3{i}=(1:m);
        end
    end
    %% 3d) Check and make sure that the upper limit of the grid is covered
    % If income (y) is large, than an agent may actually be a saver and end up
    % with more final period assets than he has at the beginning.  This isn't a
    % problem in and of itself, but we need to be careful about grids in this
    % case to ensure the upper part of the grid is covered for interpolation
    % purposes.  Here, I check if the maximum of our asset grid (max(a_grid))
    % is less than the endogenous grid we generated.  If not, then we run a
    % grid search to add in this point and ensure coverage.  I do this for each
    % of the health states.
    
    
    a_grid1=repmat(a_grid,1,hc_states); a_grid2=repmat(a_grid,1,hc_states); a_grid3=repmat(a_grid,1,hc_states);
    
    
    I1_max=(max(E_grid_1)<max(a_grid)); I2_max=(max(E_grid_2)<max(a_grid)); I3_max=(max(E_grid_3)<max(a_grid));
    
    if sum(I1_max)>0
        
        final_asset1=zeros(1,hc_states);
        maxbudget=(1+r)*max(a_grid)+y-hc1_minus;
        
        for i=1:hc_states
            a_max_grid=0:10:maxbudget(i);
            asset=max(a_grid)*ones(length(a_max_grid),1);
            [~,point] =max(Value_Function1(a_max_grid',asset,hc1_minus(i)));
            final_asset1(i)=a_max_grid(point);
            I1{i}=[I1{i};m+1];
        end
        
        E_grid_1=[E_grid_1; max(a_grid)*ones(1,hc_states)];
        a_grid1=[a_grid1; final_asset1];
    end
    
    
    if sum(I2_max)>0
        
        final_asset2=zeros(1,hc_states);
        maxbudget=(1+r)*max(a_grid)+y-hc2_minus;
        
        for i=1:hc_states
            a_max_grid=0:10:maxbudget(i);
            asset=max(a_grid)*ones(length(a_max_grid),1);
            [~,point] =max(Value_Function2(a_max_grid',asset,hc2_minus(i)));
            final_asset2(i)=a_max_grid(point);
            I2{i}=[I2{i};m+1];
        end
        E_grid_2=[E_grid_2; max(a_grid)*ones(1,hc_states)];
        a_grid2=[a_grid2; final_asset2];
    end
    
    if sum(I3_max)>0
        
        final_asset3=zeros(1,hc_states);
        maxbudget=(1+r)*max(a_grid)+y-hc3_minus;
        
        for i=1:hc_states
            a_max_grid=0:10:maxbudget(i);
            asset=max(a_grid)*ones(length(a_max_grid),1);
            [~,point] =max(Value_Function3(a_max_grid',asset,hc3_minus(i)));
            final_asset3(i)=a_max_grid(point);
            I3{i}=[I3{i};m+1];
        end
        
        E_grid_3=[E_grid_3; max(a_grid)*ones(1,hc_states)];
        a_grid3=[a_grid3; final_asset3];
    end
    %% Interpolation
    %Because as outlined in the solution algorithm there exist discrete jumps
    %in the policy functions.  To prevent interpolating over these jumps, we
    %thus only want to interpolate over those pairs that are self-insuring for
    %the same things.  As outlined in the solution algorithm, these sets
    %correspond to the connected regions of our endogenous grid, i.e. the
    %regions where we didn't throw out any points between points.  This step is
    %implemented here.
    
    %Here I save the endogenous grid pairs we kept from the previous steps.
    for i=1:hc_states
        Endog_1{i}=[E_grid_1(I1{i},i), a_grid1(I1{i},i)];
        Endog_2{i}=[E_grid_2(I2{i},i), a_grid2(I2{i},i)];
        Endog_3{i}=[E_grid_3(I3{i},i), a_grid3(I3{i},i)];
        
        %I initially interpolate over all of these pairs.
        Policy_store1{i}=interp1(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear');
        Policy_store2{i}=interp1(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear');
        Policy_store3{i}=interp1(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear');
        
    end
    
    %%
    
    
    for i=1:hc_states
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
        
               
   
        

        
        GS_I1{i}(1)=[];GS_I2{i}(1)=[];GS_I3{i}(1)=[];
        
        
        %Define the budget_bound as the maximum number of assets that one could
        %carry over to the next period.  For asset holdings below the
        %budget_bound we delete from our grid.  We do not need to solve these
        %for the policy functions since the agent has no choice but to go on
        %government care.
        
        budget_bound{i}=[find(((1+r)*a_grid+y-hc1_minus(i))>0, 1 ), ...
            find(((1+r)*a_grid+y-hc2_minus(i))>0, 1 ),...
            find(((1+r)*a_grid+y-hc3_minus(i)-chi_LTC)>0, 1 )];
        GS_I1{i}(GS_I1{i}<a_grid(budget_bound{i}(:,1)))=[];
        GS_I2{i}(GS_I2{i}<a_grid(budget_bound{i}(:,2)))=[];
        GS_I3{i}(GS_I3{i}<a_grid(budget_bound{i}(:,3)))=[];
        
        Out_of_budget1{i}=1:1:(budget_bound{i}(:,1)-1);
        Out_of_budget2{i}=1:1:(budget_bound{i}(:,2)-1);
        Out_of_budget3{i}=1:1:(budget_bound{i}(:,3)-1);
        
        if (isempty(Out_of_budget1{i})==0)
            Policy_store1{i}(Out_of_budget1{i})=0; %Set infeasible points policy to zero.
        end
        if (isempty(Out_of_budget2{i})==0)
            Policy_store2{i}(Out_of_budget2{i})=0;
        end
        if (isempty(Out_of_budget3{i})==0)
            Policy_store3{i}(Out_of_budget3{i})=0;
        end
        
        
        % Here we conduct a grid search in the regions that aren't connected
        % and store the polciy functions.  This is the big step that is a
        % correction from the Fella paper.
        if (isempty(GS_I1{i})==0)
            for asset=GS_I1{i}'
                a_next_grid=(0:.1:((1+r)*asset+y-hc1_minus(i)))';
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function1(a_next_grid,a_start,hc1_minus(i)));
                Policy_store1{i}(a_grid==asset)=a_next_grid(point);
                
            end
        end
        
        if (isempty(GS_I2{i})==0)
            for asset=GS_I2{i}'
                a_next_grid=(0:.1:((1+r)*asset+y-hc2_minus(i)))';
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function2(a_next_grid,a_start,hc2_minus(i)));
                Policy_store2{i}(a_grid==asset)=a_next_grid(point);
                
            end
        end
        
        if (isempty(GS_I3{i})==0)
            for asset=GS_I3{i}'
                a_next_grid=(0:.1:((1+r)*asset+y-hc3_minus(i)))';
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function3(a_next_grid,a_start,hc3_minus(i)));
                Policy_store3{i}(a_grid==asset)=a_next_grid(point);
            end
        end
  
        %% Now we calculate the consumption allocations implied by our policy
        % functions.
        
        Consumption_store1{i}=(1+r)*a_grid+y-hc1_minus(i)-Policy_store1{i};
        Consumption_store2{i}=(1+r)*a_grid+y-hc2_minus(i)-Policy_store2{i};
        Consumption_store3{i}=(1+r)*a_grid+y-hc3_minus(i)-Policy_store3{i};
        if (isempty(Out_of_budget1{i})==0)
            Policy_store1{i}(Out_of_budget1{i})=1e10;
        end
        if (isempty(Out_of_budget2{i})==0)
            Policy_store2{i}(Out_of_budget2{i})=1e10;
        end
        if (isempty(Out_of_budget3{i})==0)
            Policy_store3{i}(Out_of_budget3{i})=1e10;
        end
        
        
        
        
        %% Calculate Value Functions over entire grid
        Value_store1{i}=Value_Function1(Policy_store1{i}, a_grid, hc1_minus(i));
        Value_store2{i}=Value_Function2(Policy_store2{i}, a_grid, hc2_minus(i));
        Value_store3{i}=Value_Function3(Policy_store3{i}, a_grid, hc3_minus(i));
        
        
        
        
        
        %%    Check boundary condition: Check and see if the agent is better
        %     off consuming everything or sticking with the interior policy
        %     function we calculated previously.
        [Value_store1{i},Max_I1{i}]=max([Value_Function1(Policy_store1{i},a_grid, hc1_minus(i)), ...
            Value_Function1(0,a_grid, hc1_minus(i))],[],2);
        [Value_store2{i},Max_I2{i}]=max([Value_Function2(Policy_store2{i},a_grid, hc2_minus(i)), ...
            Value_Function2(0,a_grid, hc2_minus(i))],[],2);
        [Value_store3{i},Max_I3{i}]=max([Value_Function3(Policy_store3{i},a_grid, hc3_minus(i)), ...
            Value_Function3(0,a_grid, hc3_minus(i))],[],2);
        %Store the points where the agent is better off consuming
        %everything
        Boundary_I1=find(Max_I1{i}==2); Boundary_I2=find(Max_I2{i}==2); Boundary_I3=find(Max_I3{i}==2);
        
        %and set the policy function of these grid points to zero
        Policy_store1{i}(Boundary_I1)=0;
        Policy_store2{i}(Boundary_I2)=0;
        Policy_store3{i}(Boundary_I3)=0;
        
        %and the consumption policy of these grid points to all of wealth.
        Consumption_store1{i}(Boundary_I1)=(1+r)*a_grid(Boundary_I1)+y-hc1_minus(i);
        Consumption_store2{i}(Boundary_I2)=(1+r)*a_grid(Boundary_I2)+y-hc2_minus(i);
        Consumption_store3{i}(Boundary_I3)=(1+r)*a_grid(Boundary_I3)+y-hc3_minus(i);
        
        
        % Check constraint in 3rd health state
        Con_fl3{i}=find(Consumption_store3{i}<chi_LTC);
        if (isempty(Out_of_budget3{i})==0)
            Con_fl3{i}(Out_of_budget3{i}(Out_of_budget3{i}'<max(Con_fl3{i})))=[];
        end
        %If the constraint is violated then we run a grid search to see what
        %the optimal policy function is in the constrained consumption set.
        if (isempty(Con_fl3{i})==0)
            for k=Con_fl3{i}'
                asset=a_grid(k);
                a_next_grid=[(0:.1:((1+r)*asset+y-hc3_minus(i)-chi_LTC))';((1+r)*asset+y-hc3_minus(i)-chi_LTC)];
                a_start=asset'*ones(length(a_next_grid),1);
                [~,point] =max(Value_Function3(a_next_grid,a_start,hc3_minus(i)));
                Policy_store3{i}(k)=a_next_grid(point);
                Consumption_store3{i}(k)=(1+r)*asset+y-hc3_minus(i)-a_next_grid(point);
            end
        end
        
        
        %% We reassign the budget infeasible asset  and consumption policies to
        % 0 and the value function grid to be very small.  These people must go
        % on government care.
        if (isempty(Out_of_budget1{i})==0)
            Policy_store1{i}(Out_of_budget1{i})=0;
            Consumption_store1{i}(Out_of_budget1{i})=0;
            Value_store1{i}(Out_of_budget1{i})=-1e20;
        end
        if (isempty(Out_of_budget2{i})==0)
            Policy_store2{i}(Out_of_budget2{i})=0;
            Consumption_store2{i}(Out_of_budget2{i})=0;
            Value_store2{i}(Out_of_budget2{i})=-1e20;
        end
        if (isempty(Out_of_budget3{i})==0)
            Policy_store3{i}(Out_of_budget3{i})=0;
            Consumption_store3{i}(Out_of_budget3{i})=0;
            Value_store3{i}(Out_of_budget3{i})=-1e20;
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
        Gov_I1{i}=find(Value_store1{i}<gc1); Value_store1{i}(Gov_I1{i})=gc1;Policy_store1{i}(Gov_I1{i})=0; Consumption_store1{i}(Gov_I1{i})=0;
        Gov_I2{i}=find(Value_store2{i}<gc2); Value_store2{i}(Gov_I2{i})=gc2;Policy_store2{i}(Gov_I2{i})=0; Consumption_store2{i}(Gov_I2{i})=0;
        Gov_I3{i}=find(Value_store3{i}<gc3); Value_store3{i}(Gov_I3{i})=gc3;Policy_store3{i}(Gov_I3{i})=0; Consumption_store3{i}(Gov_I3{i})=0;
        
        %Here I store indexes where constraints are binding.
        Con_bind3{i}=(Consumption_store3{i}==chi_LTC);
        Gov_bind_I1{i}=(Value_store1{i}==gc1);
        Gov_bind_I2{i}=(Value_store2{i}==gc2);
        Gov_bind_I3{i}=(Value_store3{i}==gc3);
        
        %Define Value_function for state 3.  There is no decicision to be made
    %here.
        Value_store4{i}=max(bequest((1+r)*a_grid-hc4_minus(i),omega_bar,gamma,phi), bequest(0,omega_bar,gamma,phi));
        
    end
    
    %This changes the cell-arrays to matrices.  This is done so that we can
%multiply in the next steps to make it easier to calculate expectations.
    for i=1:hc_states
        Value_mat1(:,i)=Value_store1{i};
        Value_mat2(:,i)=Value_store2{i};
        Value_mat3(:,i)=Value_store3{i};
        Value_mat4(:,i)=Value_store4{i};
        Policy_mat1(:,i)=Policy_store1{i};
        Policy_mat2(:,i)=Policy_store2{i};
        Policy_mat3(:,i)=Policy_store3{i};
        Consump_mat1(:,i)=Consumption_store1{i};
        Consump_mat2(:,i)=Consumption_store2{i};
        Consump_mat3(:,i)=Consumption_store3{i};
        
        
        V_sj{1,i,T-Final_age+2-j}=Value_store1{i};
        V_sj{2,i,T-Final_age+2-j}=Value_store2{i};
        V_sj{3,i,T-Final_age+2-j}=Value_store3{i};
        V_sj{4,i,T-Final_age+2-j}=Value_store4{i};
        
        P_sj{1,i,T-Final_age+2-j}=Policy_store1{i};
        P_sj{2,i,T-Final_age+2-j}=Policy_store2{i};
        P_sj{3,i,T-Final_age+2-j}=Policy_store3{i};
      
        
        C_sj{1,i,T-Final_age+2-j}=Consumption_store1{i};
        C_sj{2,i,T-Final_age+2-j}=Consumption_store2{i};
        C_sj{3,i,T-Final_age+2-j}=Consumption_store3{i};
    end
    

end

% %% 5 - Calculate Final Values -
% %I interpolate over the value function grid and the policy function grid
% %for the specified health cost state to obtain the value function and
% %policy function for the given starting values.  I return these values as
% %output of the function.
% 
% value_mat=[interp1(a_grid, Value_mat1(:,final_h), final_a,'cubic') ,interp1(a_grid, Value_mat2(:,final_h), final_a,'cubic'),...
%     interp1(a_grid, Value_mat3(:,final_h), final_a,'cubic'), interp1(a_grid, Value_mat4(:,final_h), final_a,'cubic')];
% policy_mat=[interp1(a_grid, Policy_mat1(:,final_h), final_a,'cubic'),interp1(a_grid, Policy_mat2(:,final_h), final_a,'cubic'),...
%     interp1(a_grid, Policy_mat3(:,final_h), final_a,'cubic'), 0];
% 
% 
% V_sj=value_mat(final_s);
% P_sj=policy_mat(final_s);
% 
% 
