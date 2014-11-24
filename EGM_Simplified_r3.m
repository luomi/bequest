function [V_sj,P_sj,C_sj,CO_sj,CF_sj]   = EGM_Simplified_r3(Final_age, y_category,...
                                          T, param, gender, r, a_grid, e_states)
% and returns the Value function and Policy function (asset holdings for
% period Final_age+1) in these states.

% gender=1 if male, 0 if female
% y_category=1...5, denotes which quintile of income the person is in.
%% 1 - Initialization of parameters

sigma       = exp(param(1))+1; 
theta_LTC   = exp(param(2)); 
theta_beq   = exp(param(3));
kappa_beq   = param(4); 
psi_fam     = 1/(1+exp(param(5))); 
mu_e        = param(6); 
sigma_e     = exp(param(7));
C_F         = exp(param(8)); 
betta       = 1/(1+exp(-param(9)));  

hc_states   = 3;
% I need to make this endogenous, but for now this will work as long as it matches the healthcost.m function.

% floor_c_min=min(-chi_LTC,kappa_LTC)-1e-10*(-chi_LTC>kappa_LTC);
% 
%% Define Grids

% Fill in idiosyncratic shock grids on income and transfer threshold

kappa_grid      = nodeunif(e_states,mu_e-3.*sigma_e,mu_e+3.*sigma_e);
epsilon_prob    = normpdf(kappa_grid,mu_e,sigma_e);
epsilon_prob    = epsilon_prob/sum(epsilon_prob);

% total grid over assets and need state shifter
% sometimes we need to use individual a_grid instead
ka_grid         = gridmake(a_grid,kappa_grid); 

ma              = length(a_grid);
m               = length(ka_grid);

% Grids to be reported
V_sj        = cell(4,hc_states,T-Final_age);
P_sj        = cell(3,hc_states,T-Final_age);
C_sj        = cell(3,hc_states,T-Final_age);
CO_sj       = cell(3,hc_states,T-Final_age);
CF_sj       = cell(3,hc_states,T-Final_age);


% Preallocating all data structures I will need in this problem.
I1          = cell(1,hc_states);
I2          = cell(1,hc_states);
I3          = cell(1,hc_states);

Endog_1     = cell(1,hc_states);
Endog_2     = cell(1,hc_states);
Endog_3     = cell(1,hc_states);

Policy_store1r  = cell(1,hc_states);
Policy_store2r  = cell(1,hc_states);
Policy_store3r  = cell(1,hc_states);

Policy_store1   = zeros(m,hc_states);
Policy_store2   = zeros(m,hc_states);
Policy_store3   = zeros(m,hc_states);

I1_skip     = cell(1,hc_states);
I2_skip     = cell(1,hc_states);
I3_skip     = cell(1,hc_states);

Miss_1      = cell(1,hc_states);
Miss_2      = cell(1,hc_states);
Miss_3      = cell(1,hc_states);

lower_GS1   = cell(1,hc_states);
lower_GS2   = cell(1,hc_states);
lower_GS3   = cell(1,hc_states);

upper_GS1   = cell(1,hc_states);
upper_GS2   = cell(1,hc_states);
upper_GS3   = cell(1,hc_states);

lower_GS1A  = cell(1,hc_states);
lower_GS2A  = cell(1,hc_states);
lower_GS3A  = cell(1,hc_states);

upper_GS1A  = cell(1,hc_states);
upper_GS2A  = cell(1,hc_states);
upper_GS3A  = cell(1,hc_states);

GS_A1       = cell(1,hc_states);
GS_A2       = cell(1,hc_states);
GS_A3       = cell(1,hc_states);

budget_bound    = cell(1,hc_states);
Out_of_budget1  = cell(1,hc_states);
Out_of_budget2  = cell(1,hc_states);
Out_of_budget3  = cell(1,hc_states);

Consumption_store1  = cell(1,hc_states);
Consumption_store2  = cell(1,hc_states);
Consumption_store3  = cell(1,hc_states);

Consumptionown_store1  = cell(1,hc_states);
Consumptionown_store2  = cell(1,hc_states);
Consumptionown_store3  = cell(1,hc_states);

Consumptionfam_store1  = cell(1,hc_states);
Consumptionfam_store2  = cell(1,hc_states);
Consumptionfam_store3  = cell(1,hc_states);

Value_store1    = cell(1,hc_states);
Value_store2    = cell(1,hc_states);
Value_store3    = cell(1,hc_states);
Value_store4    = cell(1,hc_states);

Max_I1_T    = cell(1,hc_states); 
Max_I2_T    = cell(1,hc_states);
Max_I3_T    = cell(1,hc_states);

% Con_fl3     = cell(1,hc_states);

Gov_I1      = cell(1,hc_states);
Gov_I2      = cell(1,hc_states);
Gov_I3      = cell(1,hc_states);

% Con_bind3   = cell(1,hc_states);

Gov_bind_I1 = cell(1,hc_states);
Gov_bind_I2 = cell(1,hc_states);
Gov_bind_I3 = cell(1,hc_states);

EYE         = eye(4);

hc1         = zeros(hc_states,1);
hc2         = hc1;
hc3         = hc1; 
hc4         = hc1;

hc1_minus   = zeros(hc_states,1);
hc2_minus   = hc1_minus;
hc3_minus   = hc1_minus;
hc4_minus   = hc1_minus;

prob1       = zeros(hc_states,1);
prob2       = prob1;
prob3       = prob1;
prob4       = prob1;

prob1_minus = zeros(hc_states,1);
prob2_minus = prob1_minus;
prob3_minus = prob1_minus;
prob4_minus = prob1_minus;

Value_mat1      = zeros(m,hc_states);
Value_mat2      = zeros(m,hc_states);
Value_mat3      = zeros(m,hc_states);
Value_mat4      = zeros(m,hc_states);

Value_mat1_int  = zeros(ma,hc_states);
Value_mat2_int  = zeros(ma,hc_states);
Value_mat3_int  = zeros(ma,hc_states);
Value_mat4_int  = zeros(ma,hc_states);

Policy_mat1     = zeros(m,hc_states);
Policy_mat2     = zeros(m,hc_states);
Policy_mat3     = zeros(m,hc_states);

Consump_mat1    = zeros(m,hc_states);
Consump_mat2    = zeros(m,hc_states);
Consump_mat3    = zeros(m,hc_states);

deriv_1_next    = zeros(m,hc_states);
deriv_2_next    = zeros(m,hc_states);
deriv_3_next    = zeros(m,hc_states);
deriv_4_next    = zeros(m,hc_states);

cf_bind1        = cell(1,hc_states);
cf_bind2        = cell(1,hc_states);
cf_bind3        = cell(1,hc_states);

cf_bound_bind1        = cell(1,hc_states);
cf_bound_bind2        = cell(1,hc_states);
cf_bound_bind3        = cell(1,hc_states);

concave_L       = cell(1,hc_states);
concave_H       = cell(1,hc_states);

G_NC            = cell(1,hc_states);
G_C_L           = cell(1,hc_states);
G_C_H           = cell(1,hc_states);


NC_Ind          = cell(1,hc_states);
a_NC            = cell(1,hc_states);
All_Concave_I   = cell(1,hc_states);


deriv_bequestr  = zeros(ma,1);
E_grid_T1r      = zeros(ma,hc_states);
E_grid_T2r      = zeros(ma,hc_states);
E_grid_T3r      = zeros(ma,hc_states);

%% T-1 Problem
% 3a)

y   = Income(T-1,gender,y_category);

for i =1:hc_states

    % Here I generate the health costs and probabilities of each cost for time T.
    [hc1(i),prob1(i)]   = healthcost(1,T,i,gender);  
    [hc2(i),prob2(i)]   = healthcost(2,T,i,gender);
    [hc3(i),prob3(i)]   = healthcost(3,T,i,gender);
    [hc4(i),prob4(i)]   = healthcost(4,T,i,gender);
    
    % Here I generate the health costs and probabilities of each cost for time T-1.
    [hc1_minus(i),prob1_minus(i)]   = healthcost(1,T-1,i,gender); 
    [hc2_minus(i),prob2_minus(i)]   = healthcost(2,T-1,i,gender);
    [hc3_minus(i),prob3_minus(i)]   = healthcost(3,T-1,i,gender);
    [hc4_minus(i),prob4_minus(i)]   = healthcost(4,T-1,i,gender);
     
end

%% Here I define the value function for each state, 1,2,3.  It is the
% utility one gets from consumption, plus the expected value given you 
% do not die next period.  We take expectation over Health
% transition  (see multiplication by Health_Transition(T-1) matrix) and
% over iid health costs for time T (see multiplication by prob1,prob2,
% etc.).

% ltc_gridT       = squeeze(ltc_cell(:,:,T-Final_age));
ltc_gridT       = zeros(ma,hc_states);
ltc_kgridT      = repmat(ltc_gridT,e_states,1);  

% Value_Function1r         = @(c, a, h) Utility_r(c, 0, ...
%         sigma, theta_LTC)+(betta.*prob4'*bequest(max...
%         (bsxfun(@times,(1+r)*((1+r)*a'+y-c'-h),ones(hc_states,1))...
%         -bsxfun(@times,ones(1,length(a)),hc4),0),theta_beq,sigma,kappa_beq))';
% 
% Value_Function2r         = @(c, a, h) Utility_r(c, 0, ...
%         sigma, theta_LTC)+(betta.*prob4'*bequest(max...
%         (bsxfun(@times,(1+r)*((1+r)*a'+y-c'-h),ones(hc_states,1))...
%         -bsxfun(@times,ones(1,length(a)),hc4),0),theta_beq,sigma,kappa_beq))';
% 
% Value_Function3r         = @(c, a, h, ltc) Utility_r(c, 1, ...
%         sigma, theta_LTC)+(betta.*prob4'*bequest(max...
%         (bsxfun(@times,(1+r)*((1+r)*a'+y-c'-h-ltc'),ones(hc_states,1))...
%         -bsxfun(@times,ones(1,length(a)),hc4),0),theta_beq,sigma,kappa_beq))'; % ltc is of [length(a),1], or a scalar

Value_Function1         = @(co, cf, a, h, kappa) Utility(co, cf, 0, ...
        sigma, theta_LTC, kappa, psi_fam)+(betta.*prob4'*bequest(max...
        (bsxfun(@times,(1+r).*a'+y-co'-cf'-h,ones(hc_states,1))...
        -bsxfun(@times,ones(1,length(a)),hc4),0),theta_beq,sigma,kappa_beq))';

Value_Function2         = @(co, cf, a, h, kappa) Utility(co, cf, 0, ...
        sigma, theta_LTC, kappa, psi_fam)+(betta.*prob4'*bequest(max...
        (bsxfun(@times,(1+r).*a'+y-co'-cf'-h,ones(hc_states,1))...
        -bsxfun(@times,ones(1,length(a)),hc4),0),theta_beq,sigma,kappa_beq))';

Value_Function3         = @(co, cf, a, h, kappa, ltc) Utility(co, cf, 1, ...
        sigma, theta_LTC, kappa, psi_fam)+(betta.*prob4'*bequest(max...
        (bsxfun(@times,(1+r).*a'+y-co'-cf'-h-ltc',ones(hc_states,1))...
        -bsxfun(@times,ones(1,length(a)),hc4),0),theta_beq,sigma,kappa_beq))'; % ltc is of [length(a),1], or a scalar

%% T-1 Problem with Endogenous Grid
%3b) Calculate derivative of T-1 Continuation Value over asset grid.  It is
%the derivative of the bequest function if the government provided bequest
%floor is not binding, and 0 if the government bequest floor is binding.

% Note that this is the expectation already, since no health and no shifter:

deriv_bequest           = theta_beq.*(((bsxfun(@times,(1+r)*ka_grid(:,1),...
        ones(1,hc_states))-bsxfun(@times,hc4',ones(size(ka_grid(:,1)))))>0)...
        .*(1+r).*(kappa_beq+(bsxfun(@times,(1+r)*ka_grid(:,1),ones(1,hc_states))-...
        bsxfun(@times,hc4',ones(size(ka_grid(:,1)))))).^(-sigma))*prob4;

% Health_Tran_T           = Health_Transition(T-1,2-gender); %T-1 to T transition
% deriv_bequest           = deriv_bequest_vec;

% Why is deriv_bequest not smooth?  It has to do with the fact that if you
% carry positive wealth to the next state, in some states government care is are binding,
% and in other states you are not.  Since the derivative in the government
% care state is 0, one would not, as the derivative switches from 0 to
% positive (i.e., the constraint doesn't bind), then there is a jump.  This
% continues until all the jumps are finished, and then it decreases
% (decreasing marginal utility is a functtional assumption.

%% 3c) Apply endogenous grid Methodfc
% Here I generate the endogenous grid.  It comes from inverting the FOC.  I
% do this for each of the 3 states {0,1,2} where the differences emerge
% from the state dependent utility and the differences in continuation
% value.

% E_grid_T1   = (bsxfun(@times,deriv_bequest.^(-1/sigma)-y+ka_grid(:,1),ones(1,hc_states))...
%                 + bsxfun(@times,hc1_minus',ones(size(ka_grid(:,1)))))./(1+r);
% E_grid_T2   = (bsxfun(@times,deriv_bequest.^(-1/sigma)-y+ka_grid(:,1),ones(1,hc_states))...
%                 + bsxfun(@times,hc2_minus',ones(size(ka_grid(:,1)))))./(1+r);
% E_grid_T3   = (bsxfun(@times,(deriv_bequest./theta_LTC).^(-1/sigma)-y+ka_grid(:,1),ones(1,hc_states))...
%                 + bsxfun(@times,hc3_minus',ones(size(ka_grid(:,1))))+ltc_kgridT)./(1+r);

E_grid_T1   = (bsxfun(@times,(deriv_bequest./((((psi_fam^psi_fam)*((1-psi_fam)^(1-psi_fam))))^(1-sigma)))...
                .^(-1/sigma)-y+ka_grid(:,1)-ka_grid(:,2),ones(1,hc_states))...
                + bsxfun(@times,hc1_minus',ones(size(ka_grid(:,1)))))./(1+r);
E_grid_T2   = (bsxfun(@times,(deriv_bequest./((((psi_fam^psi_fam)*((1-psi_fam)^(1-psi_fam))))^(1-sigma)))...
                .^(-1/sigma)-y+ka_grid(:,1)-ka_grid(:,2),ones(1,hc_states))...
                + bsxfun(@times,hc2_minus',ones(size(ka_grid(:,1)))))./(1+r);
E_grid_T3   = (bsxfun(@times,(deriv_bequest./theta_LTC./((((psi_fam^psi_fam)*((1-psi_fam)^(1-psi_fam))))^(1-sigma)))...
                .^(-1/sigma)-y+ka_grid(:,1)-ka_grid(:,2),ones(1,hc_states))...
                + bsxfun(@times,hc3_minus',ones(size(ka_grid(:,1))))+ltc_kgridT)./(1+r);

E_grid_T1(E_grid_T1==Inf)   = 0*Inf;   
E_grid_T2(E_grid_T2==Inf)   = 0*Inf;
E_grid_T3(E_grid_T3==Inf)   = 0*Inf;

%If FOC was not invertible because the derivative of the continuation value
%was 0 then I set it to be NaN.  Note that the derivative can be zero
%because of the government care option.

%% Here I find the lower and upper bounds on the non concave regions.
% Any final grid point above concave_i_H  and anything below concave_i_L
% has a one-one mapping from the FOC, so we don't have to check for a
% global max.  Anything in G_NC_i is in a nonconcave par of the
% continuation value, and so there are multiple starting assets that map to
% the same end asset.  We thus need to check for global maxes.
%

for ee = 1:e_states
        
    E_grid_T1r      = E_grid_T1(((ee-1)*ma+1):(ee*ma),:);
    E_grid_T2r      = E_grid_T2(((ee-1)*ma+1):(ee*ma),:);
    E_grid_T3r      = E_grid_T3(((ee-1)*ma+1):(ee*ma),:);
    deriv_bequestr  = deriv_bequest(((ee-1)*ma+1):(ee*ma));
    
    for i = 1:hc_states
        [concave_L{i}, concave_H{i}]    = find_concave(deriv_bequestr);
        G_NC{i}     = (concave_L{i}:1:concave_H{i})';
        G_C_L{i}    = (1:1:(concave_L{i}-1))';
        G_C_H{i}    = ((concave_H{i}+1):1:ma)';
        NC_Ind{i}   = ismember((1:ma)',G_NC{i});
        a_NC{i}     = a_grid(G_NC{i});
        All_Concave_I{i}    = (concave_L{i}>concave_H{i});
        if (isempty(All_Concave_I{i}) == 1)
            All_Concave_I{i}    = 0;
        end
    end

    % Parallel_Indic=[(1:m*hc_states)',kron((1:hc_states)',ones(m,1))];
    NC_Indices_1    = cumsum(ismember((1:ma)',G_NC{1}));
    NC_Indices_2    = cumsum(ismember((1:ma)',G_NC{2}));
    NC_Indices_3    = cumsum(ismember((1:ma)',G_NC{3}));

    keepmat1        = zeros(ma,hc_states);
    keepmat2        = keepmat1;
    keepmat3        = keepmat1;

    for ii = 1:ma
        I1_temp     = zeros(1,hc_states);
        I2_temp     = I1_temp;
        I3_temp     = I1_temp;
        keep1       = zeros(1,hc_states);
        keep2       = keep1;
        keep3       = keep1;

        for hc_i    = 1:hc_states
            if (NC_Ind{1}(ii)==1 && All_Concave_I{1}==0)
                zi  = E_grid_T1r(ii,hc_i)*ones(length(G_NC{1}),1);
                c   = max(-kappa_grid(ee),(1+r)*zi+y-a_NC{1}-hc1_minus(hc_i));
                zi  = (y-a_NC{1}-hc1_minus(hc_i)-c)/(1+r);
                cf  = (1-psi_fam)*(c+kappa_grid(ee))-kappa_grid(ee);
                co  = c-cf;
                [~,I]           = nanmax(Value_Function1(co, cf, zi, hc1_minus(hc_i),kappa_grid(ee)));
                I1_temp(hc_i)   = I*(1-isnan(zi(1)));
                keep1(hc_i)     = (I1_temp(hc_i) == NC_Indices_1(ii));
            end

            if (NC_Ind{2}(ii)==1 && All_Concave_I{2}==0)
                zi  = E_grid_T2r(ii,hc_i)*ones(length(G_NC{2}),1);
                c   = max(-kappa_grid(ee),(1+r)*zi+y-a_NC{2}-hc2_minus(hc_i));
                zi  = (y-a_NC{2}-hc2_minus(hc_i)-c)/(1+r);
                cf  = (1-psi_fam)*(c+kappa_grid(ee))-kappa_grid(ee);
                co  = c-cf;
                [~,I]           = nanmax(Value_Function2(co, cf, zi, hc2_minus(hc_i),kappa_grid(ee)));
                I2_temp(hc_i)   = I*(1-isnan(zi(1)));
                keep2(hc_i)     = (I2_temp(hc_i)==NC_Indices_2(ii));
            end

            if (NC_Ind{3}(ii)==1 && All_Concave_I{3}==0)
                zi  = E_grid_T3r(ii,hc_i)*ones(length(G_NC{3}),1);
                c   = max(-kappa_grid(ee),(1+r)*zi+y-a_NC{3}-hc3_minus(hc_i)-ltc_gridT(ii,hc_i));
                zi  = (y-a_NC{3}-hc3_minus(hc_i)-ltc_gridT(ii,hc_i)-c)/(1+r);
                cf  = (1-psi_fam)*(c+kappa_grid(ee))-kappa_grid(ee);
                co  = c-cf;
                [~,I]           = nanmax(Value_Function3(co, cf, zi, hc3_minus(hc_i),kappa_grid(ee), ltc_gridT(ii,hc_i)));
                I3_temp(hc_i)   = I*(1-isnan(zi(1)));
                keep3(hc_i)     = (I3_temp(hc_i)==NC_Indices_3(ii));
            end

        end
        keepmat1(ii,:)  = keep1;
        keepmat2(ii,:)  = keep2;
        keepmat3(ii,:)  = keep3;
    end

    a_grid1     = repmat(a_grid,1,hc_states);
    a_grid2     = a_grid1;
    a_grid3     = a_grid1;

    for i = 1:hc_states

        if (All_Concave_I{1}==1)
            I1{i}   = (1:ma)';
        else
            I1{i}   = [G_C_L{1}; find(keepmat1(:,i));G_C_H{1}];
        end

        if (All_Concave_I{2}==1)
            I2{i}   = (1:ma)';
        else
            I2{i}   = [G_C_L{2}; find(keepmat2(:,i));G_C_H{2}];
        end
        
        if (All_Concave_I{3}==1)
            I3{i}   =(1:ma)';
        else
            I3{i}   = [G_C_L{3}; find(keepmat3(:,i));G_C_H{3}];
        end
        
        Endog_1{i}  = [E_grid_T1(I1{i},i), a_grid(I1{i})];
        Endog_2{i}  = [E_grid_T2(I2{i},i), a_grid(I2{i})];
        Endog_3{i}  = [E_grid_T3(I3{i},i), a_grid(I3{i})]; 
        

        if size(Endog_1{i}((~isnan(Endog_1{i}(:,1)) & Endog_1{i}(:,1)~=Inf & Endog_1{i}(:,1)~=-Inf),:),1)>1
            %         Policy_store1{i}=interp1_parallel(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear', 'extrap');
            Policy_store1r{i}    = interp1(Endog_1{i}((~isnan(Endog_1{i}(:,1)) & Endog_1{i}(:,1)~=Inf & Endog_1{i}(:,1)~=-Inf),1), Endog_1{i}((~isnan(Endog_1{i}(:,1)) & Endog_1{i}(:,1)~=Inf & Endog_1{i}(:,1)~=-Inf),2), a_grid,'linear', 'extrap');
        else
            Policy_store1r{i}    = a_grid.*NaN;
        end

        if size(Endog_2{i}((~isnan(Endog_2{i}(:,1)) & Endog_2{i}(:,1)~=Inf & Endog_2{i}(:,1)~=-Inf),:),1)>1
            %         Policy_store2{i}=interp1_parallel(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear', 'extrap');
            Policy_store2r{i}    = interp1(Endog_2{i}((~isnan(Endog_2{i}(:,1)) & Endog_2{i}(:,1)~=Inf & Endog_2{i}(:,1)~=-Inf),1), Endog_2{i}((~isnan(Endog_2{i}(:,1)) & Endog_2{i}(:,1)~=Inf & Endog_2{i}(:,1)~=-Inf),2), a_grid,'linear', 'extrap');
        else
            Policy_store2r{i}    = a_grid.*NaN;
        end


        if size(Endog_3{i}((~isnan(Endog_3{i}(:,1)) & Endog_3{i}(:,1)~=Inf & Endog_3{i}(:,1)~=-Inf),:),1)>1
            %         Policy_store3{i}=interp1_parallel(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear', 'extrap');
            Policy_store3r{i}    = interp1(Endog_3{i}((~isnan(Endog_3{i}(:,1)) & Endog_3{i}(:,1)~=Inf & Endog_3{i}(:,1)~=-Inf),1), Endog_3{i}((~isnan(Endog_3{i}(:,1)) & Endog_3{i}(:,1)~=Inf & Endog_3{i}(:,1)~=-Inf),2), a_grid,'linear', 'extrap');
        else
            Policy_store3r{i}    = a_grid.*NaN;
        end
        
        Policy_store1(((ee-1)*ma+1):(ee*ma),i)   = Policy_store1r{i};
        Policy_store2(((ee-1)*ma+1):(ee*ma),i)   = Policy_store2r{i};
        Policy_store3(((ee-1)*ma+1):(ee*ma),i)   = Policy_store3r{i};
    end
end

for i = 1:hc_states
    Consumption_store1{i}        = max([(1+r)*ka_grid(:,1)+y-hc1_minus(i)-Policy_store1(:,i),-ka_grid(:,2)+1],[],2);
    Consumption_store2{i}        = max([(1+r)*ka_grid(:,1)+y-hc2_minus(i)-Policy_store2(:,i),-ka_grid(:,2)+1],[],2);
    Consumption_store3{i}        = max([(1+r)*ka_grid(:,1)+y-hc3_minus(i)-Policy_store3(:,i)-ltc_kgridT(:,i),-ka_grid(:,2)+1],[],2);
    Policy_store1(:,i)           = (1+r)*ka_grid(:,1)+y-hc1_minus(i)-Consumption_store1{i};
    Policy_store2(:,i)           = (1+r)*ka_grid(:,1)+y-hc2_minus(i)-Consumption_store2{i};
    Policy_store3(:,i)           = (1+r)*ka_grid(:,1)+y-hc3_minus(i)-Consumption_store3{i}-ltc_kgridT(:,i);
end

for ee = 1:e_states
    for i = 1:hc_states

        %Here we save which members of our original grid are not included in
        %our endogenous grid.
        I1_skip{i}      = ismember(a_grid1(:,i),Endog_1{i}(:,2));
        I2_skip{i}      = ismember(a_grid2(:,i),Endog_2{i}(:,2));
        I3_skip{i}      = ismember(a_grid3(:,i),Endog_3{i}(:,2));

        %Next we store the indices of those grid points that didn't make it to
        %the endogenous grid.
        Miss_1{i}       = find(I1_skip{i}==0);
        Miss_2{i}       = find(I2_skip{i}==0);
        Miss_3{i}       = find(I3_skip{i}==0);

        %Here we take the indices of the missing points and find the lower
        %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
        %endogenous grid, but a_grid(Miss_1-1) is, then we save this point
        %here.  This will be used to define the regions we run a grid search
        %over.
        lower_GS1{i}    = ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),...
                            Endog_1{i}(:,2)).*a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i);
        lower_GS1{i}(find(-ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),Endog_1{i}(:,2))+1))    = [];


        %Here we take the indices of the missing points and find the upper
        %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
        %endogenous grid, but a_grid(Miss_1+1) is, then we save the this point
        %here.  This will be used to define the regions we run a grid search
        %over.
        upper_GS1{i}    = ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i),...
                            Endog_1{i}(:,2)).*a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i);
        upper_GS1{i}(find(-ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1))+1,i),Endog_1{i}(:,2))+1))  = [];


        %Now we save the the upper and lower bounds of the assets (not indices)
        %that we will run grid searches over.
        lower_GS1A{i}   = Endog_1{i}(ismember(Endog_1{i}(:,2),lower_GS1{i}),1);
        upper_GS1A{i}   = Endog_1{i}(ismember(Endog_1{i}(:,2),upper_GS1{i}),1);

        %The above steps are repeated for states {1,2}
        lower_GS2{i}    = ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),...
                            Endog_2{i}(:,2)).*a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i);
        lower_GS2{i}(find(-ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),Endog_2{i}(:,2))+1))    = [];

        upper_GS2{i}    = ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),...
                            Endog_2{i}(:,2)).*a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i);
        upper_GS2{i}(find(-ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),Endog_2{i}(:,2))+1))     = [];

        lower_GS2A{i}   = Endog_2{i}(ismember(Endog_2{i}(:,2),lower_GS2{i}),1);
        upper_GS2A{i}   = Endog_2{i}(ismember(Endog_2{i}(:,2),upper_GS2{i}),1);

        lower_GS3{i}    = ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),...
                            Endog_3{i}(:,2)).*a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i);
        lower_GS3{i}(find(-ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),Endog_3{i}(:,2))+1))    = [];

        upper_GS3{i}    = ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),...
                            Endog_3{i}(:,2)).*a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i);
        upper_GS3{i}(find(-ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),Endog_3{i}(:,2))+1))     = [];

        lower_GS3A{i}   = Endog_3{i}(ismember(Endog_3{i}(:,2),lower_GS3{i}),1);
        upper_GS3A{i}   = Endog_3{i}(ismember(Endog_3{i}(:,2),upper_GS3{i}),1);
        %The above steps are repeated for states {1,2}
        % This step adds in the point 0 to the lower bound asset sets if the Endogenous
        % grid we generated before stops at a positive number.  This just means
        % that we will run a grid search over [0,min_endog_grid].

        % Finally, we find all of the assets that lie within the bounds which
        % we will want to run grid searches over.

        GS_A1{i}    = 0;
        GS_A2{i}    = 0;
        GS_A3{i}    = 0;
        if (isempty(Miss_1{i})==0)
            if (isempty(lower_GS1A{i})==0)
                if (lower_GS1A{i}(end)==Endog_1{i}(end,1))
                    upper_GS1A{i}   = [upper_GS1A{i}; a_grid(end)];
                end
            end
            if (isempty(upper_GS1A{i})==0)
                if (upper_GS1A{i}(1)==Endog_1{i}(1,1))
                    lower_GS1A{i}   = [0; lower_GS1A{i}];
                end
            end

            for k = 1:min([length(lower_GS1A{i}),length(upper_GS1A{i})])
%             for k = 1:length(lower_GS1A{i})
                GS_A1{i}    = [GS_A1{i};a_grid(a_grid>=lower_GS1A{i}(k) ...
                                & a_grid<=upper_GS1A{i}(k))];
            end

        end

        if (isempty(Miss_2{i})==0)
            if (isempty(lower_GS2A{i})==0)
                if (lower_GS2A{i}(end)==Endog_2{i}(end,1))
                    upper_GS2A{i}   = [upper_GS2A{i}; a_grid(end)];
                end
            end
            if (isempty(upper_GS2A{i})==0)
                if (upper_GS2A{i}(1)==Endog_2{i}(1,1))
                    lower_GS2A{i}   = [0; lower_GS2A{i}];
                end
            end

            for k = 1:min([length(lower_GS2A{i}),length(upper_GS2A{i})])
%             for k = 1:length(lower_GS2A{i})
                GS_A2{i}    = [GS_A2{i}; a_grid(a_grid>=lower_GS2A{i}(k) ...
                                & a_grid<=upper_GS2A{i}(k))];
            end
        end

        if (isempty(Miss_3{i})==0)
            if (isempty(lower_GS3A{i})==0)
                if (lower_GS3A{i}(end)==Endog_3{i}(end,1))
                    upper_GS3A{i}   = [upper_GS3A{i}; a_grid(end)];
                end
            end
            if (isempty(upper_GS3A{i})==0)
                if (upper_GS3A{i}(1)==Endog_3{i}(1,1))
                    lower_GS3A{i}   = [0; lower_GS3A{i}];
                end
            end

            for k = 1:min([length(lower_GS3A{i}),length(upper_GS3A{i})])
%             for k = 1:length(lower_GS3A{i})
                GS_A3{i}    = [GS_A3{i}; a_grid(a_grid>=lower_GS3A{i}(k) ...
                                & a_grid<=upper_GS3A{i}(k),1)];
            end
        end



        %%

        GS_A1{i}(1)     = [];
        GS_A2{i}(1)     = [];
        GS_A3{i}(1)     = [];

        GS_A1{i}        = unique([GS_A1{i};a_grid(isnan(Policy_store1r{i}))]);
        GS_A2{i}        = unique([GS_A2{i};a_grid(isnan(Policy_store2r{i}))]);
        GS_A3{i}        = unique([GS_A3{i};a_grid(isnan(Policy_store3r{i}))]);

        %Define the budget_bound as the maximum number of assets that one could
        %carry over to the next period.  For asset holdings below the
        %budget_bound we delete from our grid.  We do not need to solve these
        %for the policy functions since the agent has no choice but to go on government care.
        budget_bound{i}    = [find(((1+r).*a_grid+y-hc1_minus(i))>0,1), ...
                              find(((1+r).*a_grid+y-hc2_minus(i))>0,1),...
                              find(((1+r).*a_grid+y-hc3_minus(i)-ltc_gridT(:,i))>0,1)];

        GS_A1{i}(GS_A1{i}<a_grid(budget_bound{i}(:,1)))     = [];
        GS_A2{i}(GS_A2{i}<a_grid(budget_bound{i}(:,2)))     = [];
        GS_A3{i}(GS_A3{i}<a_grid(budget_bound{i}(:,3)))     = [];

        Out_of_budget1{i}   = 1:1:(budget_bound{i}(:,1)-1);
        Out_of_budget2{i}   = 1:1:(budget_bound{i}(:,2)-1);
        Out_of_budget3{i}   = 1:1:(budget_bound{i}(:,3)-1);
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
    %% Conduct the grid search

    for iii = 1:ma
        asset   = a_grid(iii);
        ltc     = ltc_gridT(iii,:);
        shifter = kappa_grid(ee);
        c_lb    = max(0,-shifter);
        N_gs    = 2000;
        
        for hc_i = 1:hc_states
            if  (ismember(asset,GS_A1{hc_i})==1)
                c_gridT     = linspace(c_lb,((1+r)*asset+y-hc1_minus(hc_i)),N_gs)';
                cf_gridT    = max((1-psi_fam)*(c_gridT+shifter)-shifter,0);
                co_gridT    = c_gridT-cf_gridT;
                a_start     = bsxfun(@times,asset,ones(length(c_gridT),1));
                [~,point]   = max(Value_Function1(co_gridT,cf_gridT,a_start,hc1_minus(hc_i),shifter));
                Policy_store1r{hc_i}(iii)    = (1+r)*asset+y-hc1_minus(hc_i)-c_gridT(point);
            end


            if  (ismember(asset,GS_A2{hc_i})==1)
                c_gridT     = linspace(c_lb,((1+r)*asset+y-hc2_minus(hc_i)),N_gs)';
                cf_gridT    = max((1-psi_fam)*(c_gridT+shifter)-shifter,0);
                co_gridT    = c_gridT-cf_gridT;
                a_start     = bsxfun(@times,asset,ones(length(c_gridT),1));
                [~,point]   = max(Value_Function2(co_gridT,cf_gridT,a_start,hc2_minus(hc_i),shifter));
                Policy_store2r{hc_i}(iii)    = (1+r)*asset+y-hc2_minus(hc_i)-c_gridT(point);
            end

            if  (ismember(asset,GS_A3{hc_i})==1)
                c_gridT     = linspace(c_lb,((1+r)*asset+y-hc3_minus(hc_i)),N_gs)';
                cf_gridT    = max((1-psi_fam)*(c_gridT+shifter)-shifter,0);
                co_gridT    = c_gridT-cf_gridT;
                a_start     = bsxfun(@times,asset,ones(length(c_gridT),1));
                [~,point]   = max(Value_Function3(co_gridT,cf_gridT,a_start,hc3_minus(hc_i),shifter,ltc(hc_i)));
                Policy_store3r{hc_i}(iii)    = (1+r)*asset+y-hc3_minus(hc_i)-ltc(hc_i)-c_gridT(point);
            end
        end
    end
    
    for hc_i = 1:hc_states
        Policy_store1(((ee-1)*ma+1):(ee*ma),hc_i)   = Policy_store1r{hc_i};
        Policy_store2(((ee-1)*ma+1):(ee*ma),hc_i)   = Policy_store2r{hc_i};
        Policy_store3(((ee-1)*ma+1):(ee*ma),hc_i)   = Policy_store3r{hc_i};
    end
end

%%




for hc_i = 1:hc_states

    % Now we calculate the consumption allocations implied by our policy
    % functions.
    Consumption_store1{hc_i}        = (1+r)*ka_grid(:,1)+y-hc1_minus(hc_i)-Policy_store1(:,hc_i);
    Consumption_store2{hc_i}        = (1+r)*ka_grid(:,1)+y-hc2_minus(hc_i)-Policy_store2(:,hc_i);
    Consumption_store3{hc_i}        = (1+r)*ka_grid(:,1)+y-hc3_minus(hc_i)-Policy_store3(:,hc_i)-ltc_kgridT(:,hc_i);
    
    Consumptionown_store1{hc_i}     = psi_fam.*(Consumption_store1{hc_i} + ka_grid(:,2));
    Consumptionown_store2{hc_i}     = psi_fam.*(Consumption_store2{hc_i} + ka_grid(:,2));
    Consumptionown_store3{hc_i}     = psi_fam.*(Consumption_store3{hc_i} + ka_grid(:,2));
    
    Consumptionfam_store1{hc_i}     = (1-psi_fam).*(Consumption_store1{hc_i}+ka_grid(:,2)) - ka_grid(:,2);
    Consumptionfam_store2{hc_i}     = (1-psi_fam).*(Consumption_store2{hc_i}+ka_grid(:,2)) - ka_grid(:,2);
    Consumptionfam_store3{hc_i}     = (1-psi_fam).*(Consumption_store3{hc_i}+ka_grid(:,2)) - ka_grid(:,2);
    
    cf_bind1{hc_i}     = find(Consumptionfam_store1{hc_i}<0);
    cf_bind2{hc_i}     = find(Consumptionfam_store2{hc_i}<0);
    cf_bind3{hc_i}     = find(Consumptionfam_store3{hc_i}<0);
    
    Consumptionfam_store1{hc_i}(cf_bind1{hc_i}) = 0;
    Consumptionfam_store2{hc_i}(cf_bind2{hc_i}) = 0;
    Consumptionfam_store3{hc_i}(cf_bind3{hc_i}) = 0;
    
    Consumptionown_store1{hc_i}(cf_bind1{hc_i}) = (1+r)*ka_grid(cf_bind1{hc_i},1) ...
            +y-hc1_minus(hc_i)-Policy_store1(cf_bind1{hc_i},hc_i);
    Consumptionown_store2{hc_i}(cf_bind2{hc_i}) = (1+r)*ka_grid(cf_bind2{hc_i},1) ...
            +y-hc2_minus(hc_i)-Policy_store2(cf_bind2{hc_i},hc_i);
    Consumptionown_store3{hc_i}(cf_bind3{hc_i}) = (1+r)*ka_grid(cf_bind3{hc_i},1) ...
            +y-hc3_minus(hc_i)-Policy_store3(cf_bind3{hc_i},hc_i)-ltc_kgridT(cf_bind3{hc_i},hc_i);
        
    Consumption_store1{hc_i}(cf_bind1{hc_i}) = Consumptionown_store1{hc_i}(cf_bind1{hc_i});
    Consumption_store2{hc_i}(cf_bind2{hc_i}) = Consumptionown_store2{hc_i}(cf_bind2{hc_i});
    Consumption_store3{hc_i}(cf_bind3{hc_i}) = Consumptionown_store3{hc_i}(cf_bind3{hc_i});
    
    % Calculate Value Function values over the grids.
    Value_store1{hc_i}  = Value_Function1(Consumptionown_store1{hc_i}, Consumptionfam_store1{hc_i},ka_grid(:,1), hc1_minus(hc_i),ka_grid(:,2));
    Value_store2{hc_i}  = Value_Function2(Consumptionown_store2{hc_i}, Consumptionfam_store2{hc_i},ka_grid(:,1), hc2_minus(hc_i),ka_grid(:,2));
    Value_store3{hc_i}  = Value_Function3(Consumptionown_store3{hc_i}, Consumptionfam_store3{hc_i},ka_grid(:,1), hc3_minus(hc_i),ka_grid(:,2),ltc_kgridT(:,hc_i));
   
    
    if (isempty(Out_of_budget1{hc_i})==0)
        Consumption_store1{hc_i}(Out_of_budget1{hc_i})      = 0;
        Consumptionown_store1{hc_i}(Out_of_budget1{hc_i})   = 0;
        Consumptionfam_store1{hc_i}(Out_of_budget1{hc_i})   = 0;
    end
    if (isempty(Out_of_budget2{hc_i})==0)
        Consumption_store2{hc_i}(Out_of_budget2{hc_i})      = 0;
        Consumptionown_store2{hc_i}(Out_of_budget2{hc_i})   = 0;
        Consumptionfam_store2{hc_i}(Out_of_budget2{hc_i})   = 0;
    end
    if (isempty(Out_of_budget3{hc_i})==0)
        Consumption_store3{hc_i}(Out_of_budget3{hc_i})      = 0;
        Consumptionown_store3{hc_i}(Out_of_budget3{hc_i})   = 0;
        Consumptionfam_store3{hc_i}(Out_of_budget3{hc_i})   = 0;
    end
    
    
    

    boundary_search1    = (1:length(ka_grid(:,1)))';
    boundary_search2    = boundary_search1;
    boundary_search3    = boundary_search1;
    boundary_search1(Out_of_budget1{hc_i})  = [];
    boundary_search2(Out_of_budget2{hc_i})  = [];
    boundary_search3(Out_of_budget3{hc_i})  = [];
    
    Max_I1_T{hc_i}      = zeros(length(ka_grid(:,1)),1);
    Max_I2_T{hc_i}      = Max_I1_T{hc_i};
    Max_I3_T{hc_i}      = Max_I1_T{hc_i};
    
      
    %     Check boundary condition: Check and see if the agent is better
    %     off consuming everything or is better off with the interior policy
    %     function we calculated previously.
    [Value_store1{hc_i}(boundary_search1),Max_I1_T{hc_i}(boundary_search1)]     = ...
            max([Value_Function1(Consumption_store1{hc_i}(boundary_search1)-...
            max(0,(1-psi_fam)*(Consumption_store1{hc_i}(boundary_search1)+ka_grid(boundary_search1,2))-ka_grid(boundary_search1,2)),...
            max(0,(1-psi_fam)*(Consumption_store1{hc_i}(boundary_search1)+ka_grid(boundary_search1,2))-ka_grid(boundary_search1,2)),...
            ka_grid(boundary_search1,1),hc1_minus(hc_i),ka_grid(boundary_search1,2)),...
            Value_Function1((1+r)*ka_grid(boundary_search1,1)+y-hc1_minus(hc_i)-max(0,(1-psi_fam)*((1+r)*ka_grid(boundary_search1,1)+y-hc1_minus(hc_i)+ka_grid(boundary_search1,2))-ka_grid(boundary_search1,2)),...
            max(0,(1-psi_fam)*((1+r)*ka_grid(boundary_search1,1)+y-hc1_minus(hc_i)+ka_grid(boundary_search1,2))-ka_grid(boundary_search1,2)),...
            ka_grid(boundary_search1,1),hc1_minus(hc_i),ka_grid(boundary_search1,2))],[],2);
    [Value_store2{hc_i}(boundary_search2),Max_I2_T{hc_i}(boundary_search2)]     = ...
            max([Value_Function2(Consumption_store2{hc_i}(boundary_search2)-...
            max(0,(1-psi_fam)*(Consumption_store2{hc_i}(boundary_search2)+ka_grid(boundary_search2,2))-ka_grid(boundary_search2,2)),...
            max(0,(1-psi_fam)*(Consumption_store2{hc_i}(boundary_search2)+ka_grid(boundary_search2,2))-ka_grid(boundary_search2,2)),...
            ka_grid(boundary_search2,1),hc2_minus(hc_i),ka_grid(boundary_search2,2)),...
            Value_Function2((1+r)*ka_grid(boundary_search2,1)+y-hc2_minus(hc_i)-max(0,(1-psi_fam)*((1+r)*ka_grid(boundary_search2,1)+y-hc2_minus(hc_i)+ka_grid(boundary_search2,2))-ka_grid(boundary_search2,2)),...
            max(0,(1-psi_fam)*((1+r)*ka_grid(boundary_search2,1)+y-hc2_minus(hc_i)+ka_grid(boundary_search2,2))-ka_grid(boundary_search2,2)),...
            ka_grid(boundary_search2,1),hc2_minus(hc_i),ka_grid(boundary_search2,2))],[],2);
    [Value_store3{hc_i}(boundary_search3),Max_I3_T{hc_i}(boundary_search3)]     = ...
            max([Value_Function3(Consumption_store3{hc_i}(boundary_search3)-...
            max(0,(1-psi_fam)*(Consumption_store3{hc_i}(boundary_search3)+ka_grid(boundary_search3,2))-ka_grid(boundary_search3,2)),...
            max(0,(1-psi_fam)*(Consumption_store3{hc_i}(boundary_search3)+ka_grid(boundary_search3,2))-ka_grid(boundary_search3,2)),...
            ka_grid(boundary_search3,1),hc3_minus(hc_i),ka_grid(boundary_search3,2),ltc_kgridT(boundary_search3)),...
            Value_Function3((1+r)*ka_grid(boundary_search3,1)+y-hc3_minus(hc_i)-max(0,(1-psi_fam)*((1+r)*ka_grid(boundary_search3,1)+y-hc3_minus(hc_i)+ka_grid(boundary_search3,2))-ka_grid(boundary_search3,2)),...
            max(0,(1-psi_fam)*((1+r)*ka_grid(boundary_search3,1)+y-hc3_minus(hc_i)+ka_grid(boundary_search3,2))-ka_grid(boundary_search3,2)),...
            ka_grid(boundary_search3,1),hc3_minus(hc_i),ka_grid(boundary_search3,2),ltc_kgridT(boundary_search3))],[],2);
   
    % Store the points where the agent is better off consuming everything
    Boundary_I1_T    = find(Max_I1_T{hc_i}==2);
    Boundary_I2_T    = find(Max_I2_T{hc_i}==2);
    Boundary_I3_T    = find(Max_I3_T{hc_i}==2);
    
    %and set the policy function of these grid points to zero
    Policy_store1(Boundary_I1_T,hc_i)  = 0;
    Policy_store2(Boundary_I2_T,hc_i)  = 0;
    Policy_store3(Boundary_I3_T,hc_i)  = 0;
    
    %and the consumption policy of these grid points to all of wealth.
    Consumption_store1{hc_i}(Boundary_I1_T)  = (1+r)*ka_grid(Boundary_I1_T,1)+y-hc1_minus(hc_i);
    Consumption_store2{hc_i}(Boundary_I2_T)  = (1+r)*ka_grid(Boundary_I2_T,1)+y-hc2_minus(hc_i);
    Consumption_store3{hc_i}(Boundary_I3_T)  = (1+r)*ka_grid(Boundary_I3_T,1)+y-hc3_minus(hc_i)-ltc_kgridT(Boundary_I3_T,hc_i);
    
    Consumptionfam_store1{hc_i}(Boundary_I1_T)  = (1-psi_fam)*(Consumption_store1{hc_i}(Boundary_I1_T)+ka_grid(Boundary_I1_T,2))-ka_grid(Boundary_I1_T,2);
    Consumptionfam_store2{hc_i}(Boundary_I2_T)  = (1-psi_fam)*(Consumption_store2{hc_i}(Boundary_I2_T)+ka_grid(Boundary_I2_T,2))-ka_grid(Boundary_I2_T,2);
    Consumptionfam_store3{hc_i}(Boundary_I3_T)  = (1-psi_fam)*(Consumption_store3{hc_i}(Boundary_I3_T)+ka_grid(Boundary_I3_T,2))-ka_grid(Boundary_I3_T,2);
    
    Consumptionown_store1{hc_i}(Boundary_I1_T)  = psi_fam*(Consumption_store1{hc_i}(Boundary_I1_T)+ka_grid(Boundary_I1_T,2));
    Consumptionown_store2{hc_i}(Boundary_I2_T)  = psi_fam*(Consumption_store2{hc_i}(Boundary_I2_T)+ka_grid(Boundary_I2_T,2));
    Consumptionown_store3{hc_i}(Boundary_I3_T)  = psi_fam*(Consumption_store3{hc_i}(Boundary_I3_T)+ka_grid(Boundary_I3_T,2));
    
    cf_bound_bind1{hc_i}     = find(Consumptionfam_store1{hc_i}(Boundary_I1_T)<0);
    cf_bound_bind2{hc_i}     = find(Consumptionfam_store2{hc_i}(Boundary_I2_T)<0);
    cf_bound_bind3{hc_i}     = find(Consumptionfam_store3{hc_i}(Boundary_I3_T)<0);
    
    Consumptionfam_store1{hc_i}(cf_bound_bind1{hc_i}) = 0;
    Consumptionfam_store2{hc_i}(cf_bound_bind2{hc_i}) = 0;
    Consumptionfam_store3{hc_i}(cf_bound_bind3{hc_i}) = 0;
    
    Consumptionown_store1{hc_i}(cf_bound_bind1{hc_i}) = Consumption_store1{hc_i}(cf_bound_bind1{hc_i});
    Consumptionown_store2{hc_i}(cf_bound_bind2{hc_i}) = Consumption_store2{hc_i}(cf_bound_bind2{hc_i});
    Consumptionown_store3{hc_i}(cf_bound_bind3{hc_i}) = Consumption_store3{hc_i}(cf_bound_bind3{hc_i});
        
%     % Check constraint in 3rd health state
    
%     Con_fl3{hc_i}=find(Consumption_store3{hc_i}<-floor_c_min);
%     if (isempty(Out_of_budget3{hc_i})==0 && isempty( Con_fl3{hc_i})==0)
%         Con_fl3{hc_i}(Out_of_budget3{hc_i}<=Con_fl3{hc_i}(end))=[];
%     end
%     %If the constraint is violated then we run a grid search to see what
%     %the optimal policy function is in the constrained consumption set.
%     if (isempty(Con_fl3{hc_i})==0)
%         for jj=Con_fl3{hc_i}'
%             asset=a_grid(jj);
%             e_grid=linspace(-floor_c_min,((1+r)*asset+y-hc3_minus(hc_i)),5000)';
%             a_start=asset'*ones(length(e_grid),1);
%             [~,point] =max(Value_Function3(e_grid,a_start,hc3_minus(hc_i)));
%             Policy_store3{hc_i}(jj)=(1+r)*asset+y-hc3_minus(hc_i)-e_grid(point);
%             Consumption_store3{hc_i}(jj)=e_grid(point);
%         end
%     end
    
    
    
    % We reassign the budget infeasible asset and consumption policies to
    % 0 and the value function grid to be very small.  These people must go
    % on government care.
    if (isempty(Out_of_budget1{hc_i})==0)
        Policy_store1(Out_of_budget1{hc_i},hc_i)        = 0;
        Consumption_store1{hc_i}(Out_of_budget1{hc_i})  = 0;
        Consumptionfam_store1{hc_i}(Out_of_budget1{hc_i})  = 0;
        Consumptionown_store1{hc_i}(Out_of_budget1{hc_i})  = 0;
        Value_store1{hc_i}(Out_of_budget1{hc_i})        = -1e20;
    end
    if (isempty(Out_of_budget2{hc_i})==0)
        Policy_store2(Out_of_budget2{hc_i},hc_i)        = 0;
        Consumption_store2{hc_i}(Out_of_budget2{hc_i})  = 0;
        Consumptionfam_store3{hc_i}(Out_of_budget2{hc_i})  = 0;
        Consumptionown_store3{hc_i}(Out_of_budget2{hc_i})  = 0;
        Value_store2{hc_i}(Out_of_budget2{hc_i})        = -1e20;
    end
    if  (isempty(Out_of_budget3{hc_i})==0)
        Policy_store3(Out_of_budget3{hc_i},hc_i)        = 0;
        Consumption_store3{hc_i}(Out_of_budget3{hc_i})  = 0;
        Consumptionfam_store3{hc_i}(Out_of_budget3{hc_i})  = 0;
        Consumptionown_store3{hc_i}(Out_of_budget3{hc_i})  = 0;
        Value_store3{hc_i}(Out_of_budget3{hc_i})        = -1e20;
    end
    
    
    
    
    %% 3h) Government Care
    %Here I calculate the value of government care in each state according
    %to the allocations.
    Gov_bequest     = [prob1'*bequest(max(-hc1,0),theta_beq,sigma,kappa_beq), ...
                       prob2'*bequest(max(-hc2,0),theta_beq,sigma,kappa_beq), ...
                       prob3'*bequest(max(-hc3,0),theta_beq,sigma,kappa_beq), ...
                       prob4'*bequest(max(-hc4,0),theta_beq,sigma,kappa_beq)]';
                   
%     gc1     = Utility_r(C_F, 0, sigma,theta_LTC)+EYE(1,:)*betta*Health_Transition(T-1,2-gender)*Gov_bequest;
%     gc2     = Utility_r(C_F, 0, sigma,theta_LTC)+EYE(2,:)*betta*Health_Transition(T-1,2-gender)*Gov_bequest;
%     gc3     = Utility_r(C_F, 1, sigma,theta_LTC)+EYE(3,:)*betta*Health_Transition(T-1,2-gender)*Gov_bequest;

%     gc1     = -1e20;
%     gc2     = -1e20;
%     gc3     = -1e20;

    gc1     = @(x) max(-1e30,Utility(C_F, -min(0,x), 0, sigma,theta_LTC,x, psi_fam)+EYE(1,:)*betta*Health_Transition(T-1,2-gender)*Gov_bequest);
    gc2     = @(x) max(-1e30,Utility(C_F, -min(0,x), 0, sigma,theta_LTC,x, psi_fam)+EYE(2,:)*betta*Health_Transition(T-1,2-gender)*Gov_bequest);
    gc3     = @(x) max(-1e30,Utility(C_F, -min(0,x), 1, sigma,theta_LTC,x, psi_fam)+EYE(3,:)*betta*Health_Transition(T-1,2-gender)*Gov_bequest);
    
    %Here I Find the points on the Value function grid for which it is
    %optimal for an agent to enter government care.  These points are found
    %by comparing the value function grid to the government care option and
    %selecting the higher.  I also reassign asset and consumption policies,
    % to 0.
%     Gov_I1{hc_i}    = find(Value_store1{hc_i}<gc1);
%     Gov_I2{hc_i}    = find(Value_store2{hc_i}<gc2);
%     Gov_I3{hc_i}    = find(Value_store3{hc_i}<gc3);
%     
%     Value_store1{hc_i}(Gov_I1{hc_i})    = gc1;
%     Value_store2{hc_i}(Gov_I2{hc_i})    = gc2;
%     Value_store3{hc_i}(Gov_I3{hc_i})    = gc3;

    Gov_I1{hc_i}    = find(Value_store1{hc_i}<gc1(ka_grid(:,2)));
    Gov_I2{hc_i}    = find(Value_store2{hc_i}<gc2(ka_grid(:,2)));
    Gov_I3{hc_i}    = find(Value_store3{hc_i}<gc3(ka_grid(:,2)));
    
    Value_store1{hc_i}(Gov_I1{hc_i})    = gc1(ka_grid(Gov_I1{hc_i},2));
    Value_store2{hc_i}(Gov_I2{hc_i})    = gc2(ka_grid(Gov_I2{hc_i},2));
    Value_store3{hc_i}(Gov_I3{hc_i})    = gc3(ka_grid(Gov_I3{hc_i},2));
    
    Policy_store1(Gov_I1{hc_i},hc_i)    = 0;
    Policy_store2(Gov_I2{hc_i},hc_i)    = 0;
    Policy_store3(Gov_I3{hc_i},hc_i)    = 0;
    
    Consumption_store1{hc_i}(Gov_I1{hc_i})  = C_F;
    Consumption_store2{hc_i}(Gov_I2{hc_i})  = C_F;
    Consumption_store3{hc_i}(Gov_I3{hc_i})  = C_F;
    
    Consumptionown_store1{hc_i}(Gov_I1{hc_i})   = C_F;
    Consumptionown_store2{hc_i}(Gov_I2{hc_i})   = C_F;
    Consumptionown_store3{hc_i}(Gov_I3{hc_i})   = C_F;
    
%     Consumptionfam_store1{hc_i}(Gov_I1{hc_i})   = -min(0,ka_grid(Gov_I1{hc_i},2));
%     Consumptionfam_store2{hc_i}(Gov_I2{hc_i})   = -min(0,ka_grid(Gov_I2{hc_i},2));
%     Consumptionfam_store3{hc_i}(Gov_I3{hc_i})   = -min(0,ka_grid(Gov_I3{hc_i},2));

    Consumptionfam_store1{hc_i}(Gov_I1{hc_i})   = 0;
    Consumptionfam_store2{hc_i}(Gov_I2{hc_i})   = 0;
    Consumptionfam_store3{hc_i}(Gov_I3{hc_i})   = 0;
    
%     %Here I store indexes where constraints are binding.
%     Con_bind3{hc_i}=(Consumption_store3{hc_i}==(-floor_c_min));
%     Gov_bind_I1{hc_i}   = (Value_store1{hc_i}==gc1);
%     Gov_bind_I2{hc_i}   = (Value_store2{hc_i}==gc2);
%     Gov_bind_I3{hc_i}   = (Value_store3{hc_i}==gc3);

    Gov_bind_I1{hc_i}   = (Value_store1{hc_i}==gc1(ka_grid(:,2)));
    Gov_bind_I2{hc_i}   = (Value_store2{hc_i}==gc2(ka_grid(:,2)));
    Gov_bind_I3{hc_i}   = (Value_store3{hc_i}==gc3(ka_grid(:,2)));
    
    %Define Value_function for state 3.  There is no decicision to be made
    %here.
    Value_store4{hc_i}  = max(bequest((1+r).*ka_grid(:,1)-hc4_minus(hc_i),...
            theta_beq,sigma,kappa_beq), bequest(0,theta_beq,sigma,kappa_beq));

end

for k   = 1:hc_states
    Value_mat1(:,k)     = Value_store1{k};
    Value_mat2(:,k)     = Value_store2{k};
    Value_mat3(:,k)     = Value_store3{k};
    Value_mat4(:,k)     = Value_store4{k};
    Policy_mat1(:,k)    = Policy_store1(:,k);
    Policy_mat2(:,k)    = Policy_store2(:,k);
    Policy_mat3(:,k)    = Policy_store3(:,k);
    Consump_mat1(:,k)   = Consumption_store1{k};
    Consump_mat2(:,k)   = Consumption_store2{k};
    Consump_mat3(:,k)   = Consumption_store3{k};
    
    
    V_sj{1,k,T-Final_age}   = Value_store1{k};
    V_sj{2,k,T-Final_age}   = Value_store2{k};
    V_sj{3,k,T-Final_age}   = Value_store3{k};
    V_sj{4,k,T-Final_age}   = Value_store4{k};
    
    P_sj{1,k,T-Final_age}   = Policy_store1(:,k);
    P_sj{2,k,T-Final_age}   = Policy_store2(:,k);
    P_sj{3,k,T-Final_age}   = Policy_store3(:,k);
    
    C_sj{1,k,T-Final_age}   = Consumption_store1{k};
    C_sj{2,k,T-Final_age}   = Consumption_store2{k};
    C_sj{3,k,T-Final_age}   = Consumption_store3{k};
    
    CO_sj{1,k,T-Final_age}   = Consumptionown_store1{k};
    CO_sj{2,k,T-Final_age}   = Consumptionown_store2{k};
    CO_sj{3,k,T-Final_age}   = Consumptionown_store3{k};
    
    CF_sj{1,k,T-Final_age}   = Consumptionfam_store1{k};
    CF_sj{2,k,T-Final_age}   = Consumptionfam_store2{k};
    CF_sj{3,k,T-Final_age}   = Consumptionfam_store3{k};
end






for j   = 2:(T-Final_age)
    
    y   = Income(T-j,gender,y_category);
%     ltc_grid       = squeeze(ltc_cell(:,:,T-Final_age-j+1));
    ltc_grid       = zeros(ma,hc_states);
    ltc_kgrid      = repmat(ltc_grid,e_states,1); 

    for i = 1:hc_states
        [hc1(i),prob1(i)]   = healthcost(1,T-j+1,i,gender);% Here I generate the health costs and probabilities of each cost for time T-j+1.
        [hc2(i),prob2(i)]   = healthcost(2,T-j+1,i,gender);
        [hc3(i),prob3(i)]   = healthcost(3,T-j+1,i,gender);
        [hc4(i),prob4(i)]   = healthcost(4,T-j+1,i,gender);
        
        [hc1_minus(i),prob1_minus(i)]   = healthcost(1,T-j,i,gender); % Here I generate the health costs and probabilities of each cost for time T-j
        [hc2_minus(i),prob2_minus(i)]   = healthcost(2,T-j,i,gender);
        [hc3_minus(i),prob3_minus(i)]   = healthcost(3,T-j,i,gender);
        [hc4_minus(i),prob4_minus(i)]   = healthcost(4,T-j,i,gender);
        
        %Here I calculate the derivatives of each of the continuation values, depending on if
        %one enters state {0,1,2,3} next period.  The derivations of these expressions are presented
        %in the solution algorithm.  I enforce the relevant constraints
        %that are presented in the write-up, and each is calculated
        %analytically.
        deriv_1_next(:,i)   = (1+r)*Consumption_store1{i}.^(-sigma);
        deriv_1_next((Gov_bind_I1{i}'),i)   = 0;
        
        deriv_2_next(:,i)   = (1+r)*Consumption_store2{i}.^(-sigma);
        deriv_2_next((Gov_bind_I2{i}'),i)   = 0;
        
        deriv_3_next(:,i)   = (1+r)*theta_LTC*Consumption_store3{i}.^(-sigma);
        deriv_3_next((Gov_bind_I3{i}'),i)   = 0;
        
%         deriv_3_next((Con_bind3{i}),i)=0;
        deriv_4_next(:,i)   = max((1+r)*theta_beq.*(((kappa_beq+((1+r)*ka_grid(:,1)-hc4(i)))).^(-sigma))...
                .*(((1+r).*ka_grid(:,1)-hc4(i))>0),0);
        deriv_4_next((((1+r).*ka_grid(:,1)-hc4(i))<0),i)   = 0;
        
        
    end
    
    Health_Tran     = Health_Transition(T-j,2-gender); %Transition from T-j to T-j+1
    %Here I define the continuation values.  The function takes all
    %arguments and returns the discounted expected value of continuing in the next period
    %if you are in state {0,1,2}.  It takes expectation over health
    %transition and health cost.
    
    for k = 1:hc_states
        
        Value_mat1_int(:,k)  = reshape(Value_mat1(:,k),ma,e_states)*epsilon_prob;
        Value_mat2_int(:,k)  = reshape(Value_mat1(:,k),ma,e_states)*epsilon_prob;
        Value_mat3_int(:,k)  = reshape(Value_mat1(:,k),ma,e_states)*epsilon_prob;
        Value_mat4_int(:,k)  = reshape(Value_mat1(:,k),ma,e_states)*epsilon_prob;
        
    end
        
    Cont_1  = @(a) ECV(a, 1, a_grid,Value_mat1_int,Value_mat2_int,Value_mat3_int,Value_mat4_int,prob1,prob2,prob3,prob4, Health_Tran, betta);
    Cont_2  = @(a) ECV(a, 2, a_grid,Value_mat1_int,Value_mat2_int,Value_mat3_int,Value_mat4_int,prob1,prob2,prob3,prob4, Health_Tran, betta);
    Cont_3  = @(a) ECV(a, 3, a_grid,Value_mat1_int,Value_mat2_int,Value_mat3_int,Value_mat4_int,prob1,prob2,prob3,prob4, Health_Tran, betta);
    
    %     Defines the value function given todays asset holdings, and health
    %     cost value as a function of asset policy (a_t+1)
%     Value_Function1r     = @(a_next,a,h) Utility_r((1+r)*a+y-a_next-h,0,sigma,theta_LTC)+Cont_1(a_next)';
%     Value_Function2r     = @(a_next,a,h) Utility_r((1+r)*a+y-a_next-h,0,sigma,theta_LTC)+Cont_2(a_next)';
%     Value_Function3r     = @(a_next,a,h,ltc) Utility_r((1+r)*a+y-a_next-h-ltc,0,sigma,theta_LTC)+Cont_3(a_next)';
    
    Value_Function1     = @(a_next,a,h,kappa) Utility(max(0,psi_fam*((1+r)*a+y-a_next-h+kappa)),...
                            max(0,(1-psi_fam)*((1+r)*a+y-a_next-h+kappa)-kappa),0,sigma,theta_LTC,...
                            kappa,psi_fam) + Cont_1(a_next)';
    Value_Function2     = @(a_next,a,h,kappa) Utility(max(0,psi_fam*((1+r)*a+y-a_next-h+kappa)),...
                            max(0,(1-psi_fam)*((1+r)*a+y-a_next-h+kappa)-kappa),0,sigma,theta_LTC,...
                            kappa,psi_fam) + Cont_2(a_next)';
    Value_Function3     = @(a_next,a,h,kappa,ltc) Utility(max(0,psi_fam*((1+r)*a+y-a_next-h+kappa)),...
                            max(0,(1-psi_fam)*((1+r)*a+y-a_next-h-ltc+kappa)-kappa),0,sigma,theta_LTC,...
                            kappa,psi_fam) + Cont_3(a_next)';
    
    
    
    % Take expectation of derivative of continuation values over health
    % cost.  For next periods health state i, deriv_i_next*probi  gives the
    % expected derivative of the continuation value given you are in state
    % i.
    
    deriv_health_integrated     = [deriv_1_next*prob1, deriv_2_next*prob2, deriv_3_next*prob3, deriv_4_next*prob4];
    
    %Takes expectation over health transition.  Gives the
    % expected derivative of the continuation value given you are in state
    % i today.
    deriv_health_Cont_Val       = [Health_Tran(1,:)*deriv_health_integrated';...
                                   Health_Tran(2,:)*deriv_health_integrated';...
                                   Health_Tran(3,:)*deriv_health_integrated']';
    
%     deriv_shifter_integrated1   = reshape(deriv_health_Cont_Val(:,1),length(a_grid),e_states)*epsilon_prob;
%     deriv_shifter_integrated2   = reshape(deriv_health_Cont_Val(:,2),length(a_grid),e_states)*epsilon_prob;
%     deriv_shifter_integrated3   = reshape(deriv_health_Cont_Val(:,3),length(a_grid),e_states)*epsilon_prob;
    
    deriv_Cont_Val1             = vec(bsxfun(@times,reshape(deriv_health_Cont_Val(:,1),...
                                  length(a_grid),e_states)*epsilon_prob,ones(1,e_states)));
    deriv_Cont_Val2             = vec(bsxfun(@times,reshape(deriv_health_Cont_Val(:,2),...
                                  length(a_grid),e_states)*epsilon_prob,ones(1,e_states)));
    deriv_Cont_Val3             = vec(bsxfun(@times,reshape(deriv_health_Cont_Val(:,3),...
                                  length(a_grid),e_states)*epsilon_prob,ones(1,e_states)));
    deriv_Cont_Val              = [deriv_Cont_Val1,deriv_Cont_Val2,deriv_Cont_Val3];
    
    %% Here I generate the endogenous grid.  It comes from inverting the FOC.  I
    % do this for each of the 3 states {0,1,2} where the differences emerge
    % from the state dependent utility and the differences in continuation
    % value.
    
%     E_grid_1   = (bsxfun(@times,deriv_Cont_Val1.^(-1/sigma)-y+ka_grid(:,1),ones(1,hc_states))...
%                 + bsxfun(@times,hc1_minus',ones(size(ka_grid(:,1)))))./(1+r);
%     E_grid_2   = (bsxfun(@times,deriv_Cont_Val2.^(-1/sigma)-y+ka_grid(:,1),ones(1,hc_states))...
%                 + bsxfun(@times,hc2_minus',ones(size(ka_grid(:,1)))))./(1+r);
%     E_grid_3   = (bsxfun(@times,(deriv_Cont_Val3./theta_LTC).^(-1/sigma)-y+ka_grid(:,1),ones(1,hc_states))...
%                 + bsxfun(@times,hc3_minus',ones(size(ka_grid(:,1))))+ltc_kgrid)./(1+r);

    E_grid_1   = (bsxfun(@times,(deriv_Cont_Val1./((((psi_fam^psi_fam)*((1-psi_fam)^(1-psi_fam))))^(1-sigma)))...
                .^(-1/sigma)-y+ka_grid(:,1)-ka_grid(:,2),ones(1,hc_states))...
                + bsxfun(@times,hc1_minus',ones(size(ka_grid(:,1)))))/(1+r);
    E_grid_2   = (bsxfun(@times,(deriv_Cont_Val2./((((psi_fam^psi_fam)*((1-psi_fam)^(1-psi_fam))))^(1-sigma)))...
                .^(-1/sigma)-y+ka_grid(:,1)-ka_grid(:,2),ones(1,hc_states))...
                + bsxfun(@times,hc2_minus',ones(size(ka_grid(:,1)))))/(1+r);
    E_grid_3   = (bsxfun(@times,(deriv_Cont_Val3./theta_LTC./((((psi_fam^psi_fam)*((1-psi_fam)^(1-psi_fam))))^(1-sigma)))...
                .^(-1/sigma)-y+ka_grid(:,1)-ka_grid(:,2),ones(1,hc_states))...
                + bsxfun(@times,hc3_minus',ones(size(ka_grid(:,1))))+ltc_kgrid)/(1+r);
            
    E_grid_1(E_grid_1==Inf)     = 0*Inf;
    E_grid_2(E_grid_2==Inf)     = 0*Inf;
    E_grid_3(E_grid_3==Inf)     = 0*Inf;

    %% Here I find the lower and upper bounds on the non concave regions.
    % Any final grid point above concave_i_H  and anything below concave_i_L
    % has a one-one mapping from the FOC, so we don't have to check for a
    % global max.  Anything in G_NC_i is in a nonconcave part of the
    % continuation value, and so there are multiple starting assets that map to
    % the same end asset.  We thus need to check for global maxes.
    %
    for ee = 1:e_states
        
        E_grid_1r       = E_grid_1(((ee-1)*ma+1):(ee*ma),:);
        E_grid_2r       = E_grid_2(((ee-1)*ma+1):(ee*ma),:);
        E_grid_3r       = E_grid_3(((ee-1)*ma+1):(ee*ma),:);
        deriv_Cont_Valr = deriv_Cont_Val(((ee-1)*ma+1):(ee*ma),:);
        
        for i = 1:hc_states
            [concave_L{i}, concave_H{i}]    = find_concave(deriv_Cont_Valr(:,i));
            G_NC{i}     = (concave_L{i}:1:concave_H{i})';
            G_C_L{i}    = (1:1:(concave_L{i}-1))';
            G_C_H{i}    = ((concave_H{i}+1):1:ma)';
            NC_Ind{i}   = ismember((1:ma)',G_NC{i});
            a_NC{i}     = a_grid(G_NC{i});
            All_Concave_I{i}    = (concave_L{i}> concave_H{i});
            if (isempty(All_Concave_I{i})==1)
                All_Concave_I{i}    = 0;
            end
        end

        % Parallel_Indic=[(1:m*hc_states)',kron((1:hc_states)',ones(m,1))];
        NC_Indices_1    = cumsum(ismember((1:ma)',G_NC{1}));
        NC_Indices_2    = cumsum(ismember((1:ma)',G_NC{2}));
        NC_Indices_3    = cumsum(ismember((1:ma)',G_NC{3}));

        keepmat1    = zeros(ma,3);
        keepmat2    = keepmat1;
        keepmat3    = keepmat1;

        for iv   = 1:ma
            I1_temp     = zeros(1,hc_states);
            I2_temp     = I1_temp;
            I3_temp     = I1_temp;

            keep1       = zeros(1,hc_states);
            keep2       = keep1;
            keep3       = keep1;

            for hc_i    = 1:hc_states
            if (NC_Ind{1}(iv)==1 && All_Concave_I{1}==0)
                zi  = E_grid_1r(iv,hc_i)*ones(length(G_NC{1}),1);
                a   = a_grid(G_NC{1});
                [~,I]           = nanmax(Value_Function1(a, zi, hc1_minus(hc_i),kappa_grid(ee)));
                I1_temp(hc_i)   = I*(1-isnan(zi(1)));
                keep1(hc_i)     = (I1_temp(hc_i) == NC_Indices_1(iv));
            end

            if (NC_Ind{2}(iv)==1 && All_Concave_I{2}==0)
                zi  = E_grid_2r(iv,hc_i)*ones(length(G_NC{2}),1);
                a   = a_grid(G_NC{2});
                [~,I]           = nanmax(Value_Function2(a, zi, hc2_minus(hc_i),kappa_grid(ee)));
                I2_temp(hc_i)   = I*(1-isnan(zi(1)));
                keep2(hc_i)     = (I2_temp(hc_i)==NC_Indices_2(iv));
            end

            if (NC_Ind{3}(iv)==1 && All_Concave_I{3}==0)
                zi  = E_grid_3r(iv,hc_i)*ones(length(G_NC{3}),1);
                a   = a_grid(G_NC{3});
                [~,I]           = nanmax(Value_Function3(a, zi, hc3_minus(hc_i), kappa_grid(ee),ltc_grid(iv,hc_i)));
                I3_temp(hc_i)   = I*(1-isnan(zi(1)));
                keep3(hc_i)     = (I3_temp(hc_i)==NC_Indices_3(iv));
            end        
            keepmat1(iv,:)  = keep1;
            keepmat2(iv,:)  = keep2;
            keepmat3(iv,:)  = keep3;
            end
        end
        a_grid1     = repmat(a_grid,1,hc_states);
        a_grid2     = a_grid1;
        a_grid3     = a_grid1;

        for i = 1:hc_states

            if (All_Concave_I{1}==1)
                I1{i}   = (1:ma)';
            else
                I1{i}   = [G_C_L{1}; find(keepmat1(:,i));G_C_H{1}];
            end

            if (All_Concave_I{2}==1)
                I2{i}   = (1:ma)';
            else
                I2{i}   = [G_C_L{2}; find(keepmat2(:,i));G_C_H{2}];
            end
            if (All_Concave_I{3}==1)
                I3{i}   = (1:ma)';
            else
                I3{i}   = [G_C_L{3}; find(keepmat3(:,i));G_C_H{3}];
            end

%             Endog_1{i}  = [E_grid_1r(I1{i},i), a_grid(I1{i})];
%             Endog_2{i}  = [E_grid_2r(I2{i},i), a_grid(I2{i})];
%             Endog_3{i}  = [E_grid_3r(I3{i},i), a_grid(I3{i})]; 

            Endog_1{i}  = real([E_grid_1r(I1{i},i), a_grid(I1{i})]);
            Endog_2{i}  = real([E_grid_2r(I2{i},i), a_grid(I2{i})]);
            Endog_3{i}  = real([E_grid_3r(I3{i},i), a_grid(I3{i})]); 

            if size(Endog_1{i}((~isnan(Endog_1{i}(:,1)) & Endog_1{i}(:,1)~=Inf & Endog_1{i}(:,1)~=-Inf),:),1)>1
                %         Policy_store1{i}=interp1_parallel(Endog_1{i}(:,1) ,Endog_1{i}(:,2), a_grid,'linear', 'extrap');
                Policy_store1r{i}    = interp1(Endog_1{i}((~isnan(Endog_1{i}(:,1)) & Endog_1{i}(:,1)~=Inf & Endog_1{i}(:,1)~=-Inf),1), Endog_1{i}((~isnan(Endog_1{i}(:,1)) & Endog_1{i}(:,1)~=Inf & Endog_1{i}(:,1)~=-Inf),2), a_grid,'linear','extrap');
            else
                Policy_store1r{i}    = a_grid.*NaN;
            end

            if size(Endog_2{i}((~isnan(Endog_2{i}(:,1)) & Endog_2{i}(:,1)~=Inf & Endog_2{i}(:,1)~=-Inf),:),1)>1
                %         Policy_store2{i}=interp1_parallel(Endog_2{i}(:,1), Endog_2{i}(:,2), a_grid,'linear', 'extrap');
                Policy_store2r{i}    = interp1(Endog_2{i}((~isnan(Endog_2{i}(:,1)) & Endog_2{i}(:,1)~=Inf & Endog_2{i}(:,1)~=-Inf),1), Endog_2{i}((~isnan(Endog_2{i}(:,1)) & Endog_2{i}(:,1)~=Inf & Endog_2{i}(:,1)~=-Inf),2), a_grid,'linear', 'extrap');
            else
                Policy_store2r{i}    = a_grid.*NaN;
            end


            if size(Endog_3{i}((~isnan(Endog_3{i}(:,1)) & Endog_3{i}(:,1)~=Inf & Endog_3{i}(:,1)~=-Inf),:),1)>1
                %         Policy_store3{i}=interp1_parallel(Endog_3{i}(:,1), Endog_3{i}(:,2), a_grid,'linear', 'extrap');
                Policy_store3r{i}    = interp1(Endog_3{i}((~isnan(Endog_3{i}(:,1)) & Endog_3{i}(:,1)~=Inf & Endog_3{i}(:,1)~=-Inf),1), Endog_3{i}((~isnan(Endog_3{i}(:,1)) & Endog_3{i}(:,1)~=Inf & Endog_3{i}(:,1)~=-Inf),2), a_grid,'linear', 'extrap');
            else
                Policy_store3r{i}    = a_grid.*NaN;
            end
                
        Policy_store1(((ee-1)*ma+1):(ee*ma),i)   = Policy_store1r{i};
        Policy_store2(((ee-1)*ma+1):(ee*ma),i)   = Policy_store2r{i};
        Policy_store3(((ee-1)*ma+1):(ee*ma),i)   = Policy_store3r{i};
        end
    end


    for i = 1:hc_states
        Consumption_store1{i}        = max([(1+r)*ka_grid(:,1)+y-hc1_minus(i)-Policy_store1(:,i),-ka_grid(:,2)+1],[],2);
        Consumption_store2{i}        = max([(1+r)*ka_grid(:,1)+y-hc2_minus(i)-Policy_store2(:,i),-ka_grid(:,2)+1],[],2);
        Consumption_store3{i}        = max([(1+r)*ka_grid(:,1)+y-hc3_minus(i)-Policy_store3(:,i)-ltc_kgridT(:,i),-ka_grid(:,2)+1],[],2);
        Policy_store1(:,i)           = (1+r)*ka_grid(:,1)+y-hc1_minus(i)-Consumption_store1{i};
        Policy_store2(:,i)           = (1+r)*ka_grid(:,1)+y-hc2_minus(i)-Consumption_store2{i};
        Policy_store3(:,i)           = (1+r)*ka_grid(:,1)+y-hc3_minus(i)-Consumption_store3{i}-ltc_kgridT(:,i);
    end
    
%     for ee = 1:e_states
%         for i = 1:hc_states
%     %Here we save which members of our original grid are not included in
%             %our endogenous grid.
%             I1_skip{i}      = ismember(a_grid1(:,i),Endog_1{i}(:,2));
%             I2_skip{i}      = ismember(a_grid2(:,i),Endog_2{i}(:,2));
%             I3_skip{i}      = ismember(a_grid3(:,i),Endog_3{i}(:,2));
% 
%             %Next we store the indices of those grid points that didn't make it to
%             %the endogenous grid.
%             Miss_1{i}       = find(I1_skip{i}==0);
%             Miss_2{i}       = find(I2_skip{i}==0);
%             Miss_3{i}       = find(I3_skip{i}==0);
% 
%             %Here we take the indices of the missing points and find the lower
%             %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
%             %endogenous grid, but a_grid(Miss_1-1) is, then we save the this point
%             %here.  This will be used to define the regions we run a grid search
%             %over.
%             lower_GS1{i}    = ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),...
%                                 Endog_1{i}(:,2)).*a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i);
%             lower_GS1{i}(find(1-ismember(a_grid1(Miss_1{i}((Miss_1{i}-1)>0)-1,i),Endog_1{i}(:,2))))    =[];
% 
%             %Here we take the indices of the missing points and find the upper
%             %bound of the missing regions.  i.e., if a_grid(Miss_1) is not in the
%             %endogenous grid, but a_grid(Miss_1+1) is, then we save the this point
%             %here.  This will be used to define the regions we run a grid search
%             %over.
%             upper_GS1{i}    = ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i),...
%                                 Endog_1{i}(:,2)).*a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1(:,i)))+1,i);
%             upper_GS1{i}(find(-ismember(a_grid1(Miss_1{i}(Miss_1{i}<length(a_grid1))+1,i),Endog_1{i}(:,2))+1))  = [];
% 
%             %Now we save the the upper and lower bounds of the assets (not indices)
%             %that we will run grid searches over.
%             lower_GS1A{i}   = Endog_1{i}(ismember(Endog_1{i}(:,2),lower_GS1{i}),1);
%             upper_GS1A{i}   = Endog_1{i}(ismember(Endog_1{i}(:,2),upper_GS1{i}),1);
% 
%             %The above steps are repeated for states {1,2}
%             lower_GS2{i}    = ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),...
%                                 Endog_2{i}(:,2)).*a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i);
%             lower_GS2{i}(find(-ismember(a_grid2(Miss_2{i}((Miss_2{i}-1)>0)-1,i),Endog_2{i}(:,2))+1))    = [];
% 
%             upper_GS2{i}    = ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),...
%                                 Endog_2{i}(:,2)).*a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i);
%             upper_GS2{i}(find(-ismember(a_grid2(Miss_2{i}(Miss_2{i}<length(a_grid2(:,i)))+1,i),Endog_2{i}(:,2))+1)) = [];
% 
%             lower_GS2A{i}   = Endog_2{i}(ismember(Endog_2{i}(:,2),lower_GS2{i}),1);
%             upper_GS2A{i}   = Endog_2{i}(ismember(Endog_2{i}(:,2),upper_GS2{i}),1);
% 
%             lower_GS3{i}    = ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),...
%                                 Endog_3{i}(:,2)).*a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i);
%             lower_GS3{i}(find(-ismember(a_grid3(Miss_3{i}((Miss_3{i}-1)>0)-1,i),Endog_3{i}(:,2))+1))    = [];
% 
%             upper_GS3{i}    = ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),...
%                                 Endog_3{i}(:,2)).*a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i);
%             upper_GS3{i}(find(-ismember(a_grid3(Miss_3{i}(Miss_3{i}<length(a_grid3(:,i)))+1,i),Endog_3{i}(:,2))+1)) = [];
% 
%             lower_GS3A{i}   = Endog_3{i}(ismember(Endog_3{i}(:,2),lower_GS3{i}),1);
%             upper_GS3A{i}   = Endog_3{i}(ismember(Endog_3{i}(:,2),upper_GS3{i}),1);
%             %The above steps are repeated for states {1,2}
%             % This step adds in the point 0 to the lower bound asset sets if the Endogenous
%             % grid we generated before stops at a positive number.  This just means
%             % that we will run a grid search over [0,min_endog_grid].
% 
%             % Finally, we  find all of the assets that lie within the bounds which
%             % we will want to run grid searches over.
% 
%             GS_A1{i}    = 0;
%             GS_A2{i}    = 0;
%             GS_A3{i}    = 0;
% 
%             if (isempty(Miss_1{i})==0)
%                 if (isempty(lower_GS1A{i})==0)
%                     if (lower_GS1A{i}(end)==Endog_1{i}(end,1))
%                         upper_GS1A{i}   = [upper_GS1A{i}; a_grid(end)];
%                     end
%                 end
%                 if (isempty(upper_GS1A{i})==0)
%                     if (upper_GS1A{i}(1)==Endog_1{i}(1,1))
%                         lower_GS1A{i}   = [a_grid(1); lower_GS1A{i}];
%                     end
%                 end
% 
% %                 for k = 1:min([length(lower_GS1A{i}),length(upper_GS1A{i})])
%                 for k = 1:length(lower_GS1A{i})
%                     GS_A1{i}    = [GS_A1{i}; a_grid(a_grid>=lower_GS1A{i}(k) ...
%                                     & a_grid<=upper_GS1A{i}(k),1)];
%                 end
% 
%             end
% 
%             if (isempty(Miss_2{i})==0)
%                 if (isempty(lower_GS2A{i})==0)
%                     if (lower_GS2A{i}(end)==Endog_2{i}(end,1))
%                         upper_GS2A{i}   = [upper_GS2A{i}; a_grid(end)];
%                     end
%                 end
%                 if (isempty(upper_GS2A{i})==0)
%                     if (upper_GS2A{i}(1)==Endog_2{i}(1,1))
%                         lower_GS2A{i}   = [a_grid(1); lower_GS2A{i}];
%                     end
%                 end
% 
% %                 for k = 1:min([length(lower_GS2A{i}),length(upper_GS2A{i})])
%                 for k = 1:length(lower_GS2A{i})
%                     GS_A2{i}    = [GS_A2{i}; a_grid(a_grid>=lower_GS2A{i}(k) ...
%                                     & a_grid<=upper_GS2A{i}(k),1)];
%                 end
%             end
% 
%             if (isempty(Miss_3{i})==0)
%                 if (isempty(lower_GS3A{i})==0)
%                     if (lower_GS3A{i}(end)==Endog_3{i}(end,1))
%                         upper_GS3A{i}   = [upper_GS3A{i}; a_grid(end)];
%                     end
%                 end
%                 if (isempty(upper_GS3A{i})==0)
%                     if (upper_GS3A{i}(1)==Endog_3{i}(1,1))
%                         lower_GS3A{i}   = [a_grid(1); lower_GS3A{i}];
%                     end
%                 end
% 
% %                 for k = 1:min([length(lower_GS3A{i}),length(upper_GS3A{i})])
%                 for k = 1:length(lower_GS3A{i})
%                     GS_A3{i}    = [GS_A3{i}; a_grid(a_grid>=lower_GS3A{i}(k) ...
%                                     & a_grid<=upper_GS3A{i}(k),1)];
%                 end
%             end
% 
% 
%             %%
% 
%             GS_A1{i}(1) = [];
%             GS_A2{i}(1) = [];
%             GS_A3{i}(1) = [];
% 
%             GS_A1{i}    = unique([GS_A1{i};a_grid(isnan(Policy_store1r{i}))]);
%             GS_A2{i}    = unique([GS_A2{i};a_grid(isnan(Policy_store2r{i}))]);
%             GS_A3{i}    = unique([GS_A3{i};a_grid(isnan(Policy_store3r{i}))]);
% 
%             %Define the budget_bound as the maximum number of assets that one could
%             %carry over to the next period.  For asset holdings below the
%             %budget_bound we delete from our grid.  We do not need to solve these
%             %for the policy functions since the agent has no choice but to go on government care.
%             budget_bound{i} = [find(((1+r).*a_grid+y-hc1_minus(i))>0, 1 ), ...
%                                find(((1+r).*a_grid+y-hc2_minus(i))>0, 1 ),...
%                                find(((1+r).*a_grid+y-hc3_minus(i)-ltc_grid(:,i))>0, 1 )];
% 
%             GS_A1{i}(GS_A1{i}<a_grid(budget_bound{i}(:,1)))  = [];
%             GS_A2{i}(GS_A2{i}<a_grid(budget_bound{i}(:,2)))  = []; 
%             GS_A3{i}(GS_A3{i}<a_grid(budget_bound{i}(:,3)))  = [];
% 
%             Out_of_budget1{i}   =1:1:(budget_bound{i}(:,1)-1);
%             Out_of_budget2{i}   =1:1:(budget_bound{i}(:,2)-1);
%             Out_of_budget3{i}   =1:1:(budget_bound{i}(:,3)-1);
% 
%         %     if ((Endog_1{i}(end,1)<a_grid(end)) && ismember(Endog_1{i}(end,1),GS_A1{i})==0)
%         %         GS_A1{i}=[GS_A1{i},a_grid(end)];
%         %     end
%         %
%         %     if ((Endog_2{i}(end,1)<a_grid(end)) && ismember(Endog_2{i}(end,1),GS_A2{i})==0)
%         %         GS_A2{i}=[GS_A2{i},a_grid(end)];
%         %     end
%         %     if ((Endog_3{i}(end,1)<a_grid(end)) && ismember(Endog_3{i}(end,1),GS_A3{i})==0)
%         %         GS_A3{i}=[GS_A3{i},a_grid(end)];
%         %     end
%         end
%         %%
% 
%         for v = 1:ma
%             asset   = a_grid(v);
%             N_gs    = 2000;
%             ltc     = ltc_gridT(v,:);
%             
%             for hc_i = 1:hc_states
%                 if  (ismember(asset,GS_A1{hc_i})==1)
%                     a_next_grid = linspace(0,(1+r)*asset+y-hc1_minus(hc_i)-max(0,-kappa_grid(ee)),N_gs)';
%                     a_start     = bsxfun(@times,asset,ones(length(a_next_grid),1));
%                     [~,point]   = max(Value_Function1(a_next_grid,a_start,hc1_minus(hc_i),kappa_grid(ee)));
%                     Policy_store1r{hc_i}(v)    = a_next_grid(point);
%                 end
% 
% 
%                 if  (ismember(asset,GS_A2{hc_i})==1)
%                     a_next_grid = linspace(0,(1+r)*asset+y-hc2_minus(hc_i)-max(0,-kappa_grid(ee)),N_gs)';
%                     a_start     = bsxfun(@times,asset,ones(length(a_next_grid),1));
%                     [~,point]   = max(Value_Function2(a_next_grid,a_start,hc2_minus(hc_i),kappa_grid(ee)));
%                     Policy_store2r{hc_i}(v)    = a_next_grid(point);
%                 end
% 
%                 if  (ismember(asset,GS_A3{hc_i})==1)
%                     a_next_grid = linspace(0,(1+r)*asset+y-hc3_minus(hc_i)-ltc(hc_i)-max(0,-kappa_grid(ee)),N_gs)';
%                     a_start     = bsxfun(@times,asset,ones(length(a_next_grid),1));
%                     [~,point]   = max(Value_Function3(a_next_grid,a_start,hc3_minus(hc_i),kappa_grid(ee),ltc(hc_i)));
%                     Policy_store3r{hc_i}(v)    = a_next_grid(point);
%                 end
%             end
%         end
%         
%         for hc_i = 1:hc_states
%             Policy_store1(((ee-1)*ma+1):(ee*ma),hc_i)   = Policy_store1r{hc_i};
%             Policy_store2(((ee-1)*ma+1):(ee*ma),hc_i)   = Policy_store2r{hc_i};
%             Policy_store3(((ee-1)*ma+1):(ee*ma),hc_i)   = Policy_store3r{hc_i};
%         end
%     end

    
    %%
    for hc_i=1:hc_states      
        %% Now we calculate the consumption allocations implied by our policy
        % functions.
        Consumption_store1{hc_i}        = (1+r)*ka_grid(:,1)+y-hc1_minus(hc_i)-Policy_store1(:,hc_i);
        Consumption_store2{hc_i}        = (1+r)*ka_grid(:,1)+y-hc2_minus(hc_i)-Policy_store2(:,hc_i);
        Consumption_store3{hc_i}        = (1+r)*ka_grid(:,1)+y-hc3_minus(hc_i)-ltc_kgrid(:,hc_i)-Policy_store3(:,hc_i);

        Consumptionown_store1{hc_i}     = psi_fam.*(Consumption_store1{hc_i} + ka_grid(:,2));
        Consumptionown_store2{hc_i}     = psi_fam.*(Consumption_store2{hc_i} + ka_grid(:,2));
        Consumptionown_store3{hc_i}     = psi_fam.*(Consumption_store3{hc_i} + ka_grid(:,2));

        Consumptionfam_store1{hc_i}     = (1-psi_fam).*(Consumption_store1{hc_i}+ka_grid(:,2)) - ka_grid(:,2);
        Consumptionfam_store2{hc_i}     = (1-psi_fam).*(Consumption_store2{hc_i}+ka_grid(:,2)) - ka_grid(:,2);
        Consumptionfam_store3{hc_i}     = (1-psi_fam).*(Consumption_store3{hc_i}+ka_grid(:,2)) - ka_grid(:,2);

        cf_bind1{hc_i}     = find(Consumptionfam_store1{hc_i}<0);
        cf_bind2{hc_i}     = find(Consumptionfam_store2{hc_i}<0);
        cf_bind3{hc_i}     = find(Consumptionfam_store3{hc_i}<0);

        Consumptionfam_store1{hc_i}(cf_bind1{hc_i}) = 0;
        Consumptionfam_store2{hc_i}(cf_bind2{hc_i}) = 0;
        Consumptionfam_store3{hc_i}(cf_bind3{hc_i}) = 0;

        Consumptionown_store1{hc_i}(cf_bind1{hc_i}) = Consumptionown_store1{hc_i}(cf_bind1{hc_i});
        Consumptionown_store2{hc_i}(cf_bind2{hc_i}) = Consumptionown_store2{hc_i}(cf_bind2{hc_i});
        Consumptionown_store3{hc_i}(cf_bind3{hc_i}) = Consumptionown_store3{hc_i}(cf_bind3{hc_i});
        
        % Calculate Value Function values over the grids.
        Value_store1{hc_i}  = Value_Function1(Policy_store1(:,hc_i), ka_grid(:,1), hc1_minus(hc_i),ka_grid(:,2));
        Value_store2{hc_i}  = Value_Function2(Policy_store2(:,hc_i), ka_grid(:,1), hc2_minus(hc_i),ka_grid(:,2));
        Value_store3{hc_i}  = Value_Function3(Policy_store3(:,hc_i), ka_grid(:,1), hc3_minus(hc_i),ka_grid(:,2),ltc_kgrid(:,hc_i));

        
        if (isempty(Out_of_budget1{hc_i})==0)
            Consumption_store1{hc_i}(Out_of_budget1{hc_i})=0;
            Consumptionown_store1{hc_i}(Out_of_budget1{hc_i})=0;
            Consumptionfam_store1{hc_i}(Out_of_budget1{hc_i})=0;
        end
        if (isempty(Out_of_budget2{hc_i})==0)
            Consumption_store2{hc_i}(Out_of_budget2{hc_i})=0;
            Consumptionown_store2{hc_i}(Out_of_budget2{hc_i})=0;
            Consumptionfam_store2{hc_i}(Out_of_budget2{hc_i})=0;
        end
        if (isempty(Out_of_budget3{hc_i})==0)
            Consumption_store3{hc_i}(Out_of_budget3{hc_i})=0;
            Consumptionown_store3{hc_i}(Out_of_budget3{hc_i})=0;
            Consumptionfam_store3{hc_i}(Out_of_budget3{hc_i})=0;
        end
        
        
        

            
%         boundary_search1    = (1:length(ka_grid(:,1)))';
%         boundary_search2    = boundary_search1;
%         boundary_search3    = boundary_search1;
%         boundary_search1(Out_of_budget1{hc_i})  = [];
%         boundary_search2(Out_of_budget2{hc_i})  = [];
%         boundary_search3(Out_of_budget3{hc_i})  = [];
% 
%         Max_I1_T{hc_i}      = zeros(length(ka_grid(:,1)),1);
%         Max_I2_T{hc_i}      = Max_I1_T{hc_i};
%         Max_I3_T{hc_i}      = Max_I1_T{hc_i};
%     
%       
%         
%         
%         
%         %     Check boundary condition: Check and see if the agent is better
%         %     off consuming everything or s king with the interior policy
%         %     function we calculated previously.
%         
%         [Value_store1{hc_i}(boundary_search1),Max_I1_T{hc_i}(boundary_search1)]     = ...
%             max([Value_Function1(Policy_store1(boundary_search1,hc_i),...
%             ka_grid(boundary_search1,1),hc1_minus(hc_i),ka_grid(boundary_search1,2)), ...
%             Value_Function1(0,ka_grid(boundary_search1,1), hc1_minus(hc_i),ka_grid(boundary_search1,2))],[],2);
% 
%         [Value_store2{hc_i}(boundary_search2),Max_I2_T{hc_i}(boundary_search2)]     = ...
%             max([Value_Function2(Policy_store2(boundary_search2,hc_i),...
%             ka_grid(boundary_search2,1),hc2_minus(hc_i),ka_grid(boundary_search2,2)), ...
%             Value_Function2(0,ka_grid(boundary_search2,1), hc2_minus(hc_i),ka_grid(boundary_search2,2))],[],2);
% 
%         [Value_store3{hc_i}(boundary_search3),Max_I3_T{hc_i}(boundary_search3)]     = ...
%             max([Value_Function3(Policy_store3(boundary_search3,hc_i),...
%             ka_grid(boundary_search3,1),hc3_minus(hc_i), ka_grid(boundary_search3,2),ltc_kgrid(boundary_search3,hc_i)),...
%             Value_Function3(0,ka_grid(boundary_search3,1), hc3_minus(hc_i),ka_grid(boundary_search3,2), ...
%             ltc_kgrid(boundary_search3,hc_i))],[],2);
% 
%         
%         % Store the points where the agent is better off consuming everything
%         Boundary_I1_T    = find(Max_I1_T{hc_i}==2);
%         Boundary_I2_T    = find(Max_I2_T{hc_i}==2);
%         Boundary_I3_T    = find(Max_I3_T{hc_i}==2);
% 
%         %and set the policy function of these grid points to zero
%         Policy_store1(Boundary_I1_T,hc_i)  = 0;
%         Policy_store2(Boundary_I2_T,hc_i)  = 0;
%         Policy_store3(Boundary_I3_T,hc_i)  = 0;
% 
%         %and the consumption policy of these grid points to all of wealth.
%         Consumption_store1{hc_i}(Boundary_I1_T)  = (1+r)*ka_grid(Boundary_I1_T,1)+y-hc1_minus(hc_i);
%         Consumption_store2{hc_i}(Boundary_I2_T)  = (1+r)*ka_grid(Boundary_I2_T,1)+y-hc2_minus(hc_i);
%         Consumption_store3{hc_i}(Boundary_I3_T)  = (1+r)*ka_grid(Boundary_I3_T,1)+y-hc3_minus(hc_i)-ltc_kgrid(Boundary_I3_T,hc_i);
% 
%         Consumptionfam_store1{hc_i}(Boundary_I1_T)  = (1-psi_fam)*(Consumption_store1{hc_i}(Boundary_I1_T)+ka_grid(Boundary_I1_T,2))-ka_grid(Boundary_I1_T,2);
%         Consumptionfam_store2{hc_i}(Boundary_I2_T)  = (1-psi_fam)*(Consumption_store2{hc_i}(Boundary_I2_T)+ka_grid(Boundary_I2_T,2))-ka_grid(Boundary_I2_T,2);
%         Consumptionfam_store3{hc_i}(Boundary_I3_T)  = (1-psi_fam)*(Consumption_store3{hc_i}(Boundary_I3_T)+ka_grid(Boundary_I3_T,2))-ka_grid(Boundary_I3_T,2);
% 
%         Consumptionown_store1{hc_i}(Boundary_I1_T)  = psi_fam*(Consumption_store1{hc_i}(Boundary_I1_T)+ka_grid(Boundary_I1_T,2));
%         Consumptionown_store2{hc_i}(Boundary_I2_T)  = psi_fam*(Consumption_store2{hc_i}(Boundary_I2_T)+ka_grid(Boundary_I2_T,2));
%         Consumptionown_store3{hc_i}(Boundary_I3_T)  = psi_fam*(Consumption_store3{hc_i}(Boundary_I3_T)+ka_grid(Boundary_I3_T,2));
% 
%         cf_bound_bind1{hc_i}     = find(Consumptionfam_store1{hc_i}(Boundary_I1_T)<0);
%         cf_bound_bind2{hc_i}     = find(Consumptionfam_store2{hc_i}(Boundary_I2_T)<0);
%         cf_bound_bind3{hc_i}     = find(Consumptionfam_store3{hc_i}(Boundary_I3_T)<0);
%     
%         Consumptionfam_store1{hc_i}(cf_bound_bind1{hc_i}) = 0;
%         Consumptionfam_store2{hc_i}(cf_bound_bind2{hc_i}) = 0;
%         Consumptionfam_store3{hc_i}(cf_bound_bind3{hc_i}) = 0;
% 
%         Consumptionown_store1{hc_i}(cf_bound_bind1{hc_i}) = Consumption_store1{hc_i}(cf_bound_bind1{hc_i});
%         Consumptionown_store2{hc_i}(cf_bound_bind2{hc_i}) = Consumption_store2{hc_i}(cf_bound_bind2{hc_i});
%         Consumptionown_store3{hc_i}(cf_bound_bind3{hc_i}) = Consumption_store3{hc_i}(cf_bound_bind3{hc_i});

%         % Check constraint in 3rd health state
%         
%         Con_fl3{hc_i}=find(Consumption_store3{hc_i}<-floor_c_min);
%         if (isempty(Out_of_budget3{hc_i})==0  && isempty( Con_fl3{hc_i})==0)
%             Con_fl3{hc_i}(Out_of_budget3{hc_i}<=max(Con_fl3{hc_i}))=[];
%         end
%         %If the constraint is violated then we run a grid search to see what
%         %the optimal policy function is in the constrained consumption set.
%         if (isempty(Con_fl3{hc_i})==0)
%             for k=Con_fl3{hc_i}'
%                 asset=a_grid(k);
%                 a_next_grid=linspace(0,((1+r)*asset+y-hc3_minus(i)+floor_c_min),5000)';
%                 a_start=asset'*ones(length(a_next_grid),1);
%                 [~,point] =max(Value_Function3(a_next_grid,a_start,hc3_minus(hc_i)));
%                 Policy_store3{hc_i}(k)=a_next_grid(point);
%                 Consumption_store3{hc_i}(k)=(1+r)*asset+y-hc3_minus(hc_i)-a_next_grid(point);
%             end
%         end
        
        
        
        % We reassign the budget infeasible asset  and consumption policies to
        % 0 and the value function grid to be very small.  These people must go
        % on government care.
        if (isempty(Out_of_budget1{hc_i})==0)
            Policy_store1(Out_of_budget1{hc_i},hc_i)        = 0;
            Consumption_store1{hc_i}(Out_of_budget1{hc_i})  = 0;
            Consumptionown_store1{hc_i}(Out_of_budget1{hc_i})   = 0;
            Consumptionfam_store1{hc_i}(Out_of_budget1{hc_i})   = 0;
            Value_store1{hc_i}(Out_of_budget1{hc_i})        = -1e20;
        end
        if (isempty(Out_of_budget2{hc_i})==0)
            Policy_store2(Out_of_budget2{hc_i},hc_i)        = 0;
            Consumption_store2{hc_i}(Out_of_budget2{hc_i})  = 0;
            Consumptionown_store2{hc_i}(Out_of_budget2{hc_i})   = 0;
            Consumptionfam_store2{hc_i}(Out_of_budget2{hc_i})   = 0;
            Value_store2{hc_i}(Out_of_budget2{hc_i})        = -1e20;
        end
        if  (isempty(Out_of_budget3{hc_i})==0)
            Policy_store3(Out_of_budget3{hc_i},hc_i)        = 0;
            Consumption_store3{hc_i}(Out_of_budget3{hc_i})  = 0;
            Consumptionown_store3{hc_i}(Out_of_budget3{hc_i})   = 0;
            Consumptionfam_store3{hc_i}(Out_of_budget3{hc_i})   = 0;
            Value_store3{hc_i}(Out_of_budget3{hc_i})        = -1e20;
        end
        
        
        
        
        %% 3h) Government Care
        %Here I calculate the value of government care in each state according
        %to the allocations.
        
%         gc1     = Utility_r(C_F, 0,sigma,theta_LTC)+Cont_1(0)';
%         gc2     = Utility_r(C_F, 0,sigma,theta_LTC)+Cont_2(0)';
%         gc3     = Utility_r(C_F, 1,sigma,theta_LTC)+Cont_3(0)';

%         gc1     = -1e20;
%         gc2     = -1e20;
%         gc3     = -1e20;
        
        gc1     = @(x) max(-1e100,Utility(C_F, -min(0,x),0,sigma,theta_LTC,x, psi_fam)+Cont_1(0)');
        gc2     = @(x) max(-1e100,Utility(C_F, -min(0,x),0,sigma,theta_LTC,x, psi_fam)+Cont_2(0)');
        gc3     = @(x) max(-1e100,Utility(C_F, -min(0,x),1,sigma,theta_LTC,x, psi_fam)+Cont_3(0)');        
        
        %Here I Find the points on the Value function grid for which it is
        %optimal for an agent to enter government care.  These points are found
        %by comparing the value function grid to the government care option and
        %selecting the higher.  I also reassign asset and consumption policies,
        % to 0.

%         Gov_I1{hc_i}    = find(Value_store1{hc_i}<gc1);
%         Gov_I2{hc_i}    = find(Value_store2{hc_i}<gc2);
%         Gov_I3{hc_i}    = find(Value_store3{hc_i}<gc3);
% 
%         Value_store1{hc_i}(Gov_I1{hc_i})    = gc1;
%         Value_store2{hc_i}(Gov_I2{hc_i})    = gc2;
%         Value_store3{hc_i}(Gov_I3{hc_i})    = gc3;

        Gov_I1{hc_i}    = find(Value_store1{hc_i}<gc1(ka_grid(:,2)));
        Gov_I2{hc_i}    = find(Value_store2{hc_i}<gc2(ka_grid(:,2)));
        Gov_I3{hc_i}    = find(Value_store3{hc_i}<gc3(ka_grid(:,2)));

        Value_store1{hc_i}(Gov_I1{hc_i})    = gc1(ka_grid(Gov_I1{hc_i},2));
        Value_store2{hc_i}(Gov_I2{hc_i})    = gc2(ka_grid(Gov_I2{hc_i},2));
        Value_store3{hc_i}(Gov_I3{hc_i})    = gc3(ka_grid(Gov_I3{hc_i},2));
        
        Policy_store1(Gov_I1{hc_i},hc_i)    = 0;
        Policy_store2(Gov_I2{hc_i},hc_i)    = 0;
        Policy_store3(Gov_I3{hc_i},hc_i)    = 0;

        Consumption_store1{hc_i}(Gov_I1{hc_i})  = C_F;
        Consumption_store2{hc_i}(Gov_I2{hc_i})  = C_F;
        Consumption_store3{hc_i}(Gov_I3{hc_i})  = C_F;

        Consumptionown_store1{hc_i}(Gov_I1{hc_i})   = C_F;
        Consumptionown_store2{hc_i}(Gov_I2{hc_i})   = C_F;
        Consumptionown_store3{hc_i}(Gov_I3{hc_i})   = C_F;

%         Consumptionfam_store1{hc_i}(Gov_I1{hc_i})   = -min(0,ka_grid(Gov_I1{hc_i},2));
%         Consumptionfam_store2{hc_i}(Gov_I2{hc_i})   = -min(0,ka_grid(Gov_I2{hc_i},2));
%         Consumptionfam_store3{hc_i}(Gov_I3{hc_i})   = -min(0,ka_grid(Gov_I3{hc_i},2));

        Consumptionfam_store1{hc_i}(Gov_I1{hc_i})   = 0;
        Consumptionfam_store2{hc_i}(Gov_I2{hc_i})   = 0;
        Consumptionfam_store3{hc_i}(Gov_I3{hc_i})   = 0;
        
        %Here I store indexes where constraints are binding.
%         Con_bind3{hc_i}=(Consumption_store3{hc_i}==chi_LTC);
%         Gov_bind_I1{hc_i}   = (Value_store1{hc_i}==gc1);
%         Gov_bind_I2{hc_i}   = (Value_store2{hc_i}==gc2);
%         Gov_bind_I3{hc_i}   = (Value_store3{hc_i}==gc3);

        Gov_bind_I1{hc_i}   = (Value_store1{hc_i}==gc1(ka_grid(:,2)));
        Gov_bind_I2{hc_i}   = (Value_store2{hc_i}==gc2(ka_grid(:,2)));
        Gov_bind_I3{hc_i}   = (Value_store3{hc_i}==gc3(ka_grid(:,2)));
        %Define Value_function for state 3.  There is no decicision to be made
        %here.
        Value_store4{hc_i}  = max(bequest((1+r).*ka_grid(:,1)-hc4_minus(hc_i),...
                theta_beq,sigma,kappa_beq), bequest(0,theta_beq,sigma,kappa_beq));

        
            
            
            
      
        
    end
    
    for k   = 1:hc_states
        
        Value_mat1(:,k)     = Value_store1{k};
        Value_mat2(:,k)     = Value_store2{k};
        Value_mat3(:,k)     = Value_store3{k};
        Value_mat4(:,k)     = Value_store4{k};
        Policy_mat1(:,k)    = Policy_store1(:,k);
        Policy_mat2(:,k)    = Policy_store2(:,k);
        Policy_mat3(:,k)    = Policy_store3(:,k);
        Consump_mat1(:,k)   = Consumption_store1{k};
        Consump_mat2(:,k)   = Consumption_store2{k};
        Consump_mat3(:,k)   = Consumption_store3{k};
        
        
        V_sj{1,k,T-Final_age+1-j}   = Value_store1{k};
        V_sj{2,k,T-Final_age+1-j}   = Value_store2{k};
        V_sj{3,k,T-Final_age+1-j}   = Value_store3{k};
        V_sj{4,k,T-Final_age+1-j}   = Value_store4{k};
        
        P_sj{1,k,T-Final_age+1-j}   = Policy_store1(:,k);
        P_sj{2,k,T-Final_age+1-j}   = Policy_store2(:,k);
        P_sj{3,k,T-Final_age+1-j}   = Policy_store3(:,k);
        
        
        C_sj{1,k,T-Final_age+1-j}   = Consumption_store1{k};
        C_sj{2,k,T-Final_age+1-j}   = Consumption_store2{k};
        C_sj{3,k,T-Final_age+1-j}   = Consumption_store3{k};
        
        CO_sj{1,k,T-Final_age+1-j}  = Consumptionown_store1{k};
        CO_sj{2,k,T-Final_age+1-j}  = Consumptionown_store2{k};
        CO_sj{3,k,T-Final_age+1-j}  = Consumptionown_store3{k};

        CF_sj{1,k,T-Final_age+1-j}  = Consumptionfam_store1{k};
        CF_sj{2,k,T-Final_age+1-j}  = Consumptionfam_store2{k};
        CF_sj{3,k,T-Final_age+1-j}  = Consumptionfam_store3{k};
    end
   

end

end