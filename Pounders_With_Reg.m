function out=Pounders_With_Reg(in_func, initialpoint, delta_max, eta_0,eta_1, gamma_0, gamma_1, delta_0, lb, ub,weights, ParallelNodes)
%%  Initialize values
options = optimset('TolCon',1e-12,'TolFun',1e-12,'TolX',1e-12 ,'MaxFunEvals',1e5 ,'MaxIter', 1e5);
Param_N=length(initialpoint);
delta=delta_0;
x_store=initialpoint';
xk_store=initialpoint';
counter=1;
norm=@(x,y) max(weights.*(x-y));
EvalPoints=zeros(ParallelNodes,1);
%% Initial Problem
% matlabpool(ParallelNodes)
disp(1)
grid_lb=max(xk_store-delta*weights.^(-1),lb);
grid_ub=min(xk_store+delta*weights.^(-1),ub);
range=grid_ub-grid_lb;
gridpoints=(repmat(range, ParallelNodes,1).*rand(ParallelNodes, Param_N)+repmat(grid_lb,ParallelNodes,1));
parfor i=1:ParallelNodes;
    if i==1
        EvalPoints(i)=in_func(xk_store);
    end
    if i>1
        EvalPoints(i)=in_func(gridpoints(i,:));
    end
end
gridpoints(1,:)=[];
f_of_xk_store=EvalPoints(1);
model_reg=QuadraticRegression([xk_store; gridpoints],EvalPoints);
Constant=model_reg.Constant;
Linear=model_reg.Linear;
Squared=model_reg.Squared;
polynomial_k=@(parameters) Constant+parameters*Linear+(1/2)*parameters*Squared*parameters';
xk_plus_sk_store=fmincon(polynomial_k, xk_store, [], [] ,[], [], grid_lb, grid_ub,[], options);

%%

loopcount=0;
tol=1e3;
while (loopcount<250 && tol>1e-9)
    loopcount=loopcount+1;
    grid_lb=max(xk_plus_sk_store-delta*weights.^(-1),lb);
    grid_ub=min(xk_plus_sk_store+delta*weights.^(-1),ub);
    range=grid_ub-grid_lb;
    
    skippoints=1;
    gridpoints=(repmat(range, ParallelNodes,1).*rand(ParallelNodes, Param_N)+repmat(grid_lb,ParallelNodes,1));
    
    EvalPoints=zeros(ParallelNodes,1);
    parfor i=1:ParallelNodes;
        
        if i==1;
            EvalPoints(i)=in_func(xk_plus_sk_store);
        elseif i>skippoints
            EvalPoints(i)=in_func(gridpoints(i,:));
        end
    end
    gridpoints(1:skippoints,:)=[];
    f_of_xk_plus_sk_store=EvalPoints(1);
    
    rho_k=(f_of_xk_store-f_of_xk_plus_sk_store)/(polynomial_k(xk_store)-polynomial_k(xk_plus_sk_store));
    
    if rho_k<eta_0
        delta=gamma_0*delta;
        
        
        grid_lb=max(xk_store-delta*weights.^(-1),lb);
        grid_ub=min(xk_store+delta*weights.^(-1),ub);
        range=grid_ub-grid_lb;
        gridpoints=(repmat(range, ParallelNodes-1,1).*rand(ParallelNodes-1, Param_N)+repmat(grid_lb,ParallelNodes-1,1));
        EvalPoints(1)=f_of_xk_store;
        parfor i=1:(ParallelNodes-1);
            EvalPoints(i+1)=in_func(gridpoints(i,:));
        end
        model_reg=QuadraticRegression([xk_store; gridpoints],EvalPoints);
        Constant=model_reg.Constant;
        Linear=model_reg.Linear;
        Squared=model_reg.Squared;
        polynomial_k=@(parameters) Constant+parameters*Linear+(1/2)*parameters*Squared*parameters';
        xk_plus_sk_store=fmincon(polynomial_k, xk_store, [], [] ,[], [], grid_lb, grid_ub,[], options);
        
    elseif rho_k>=eta_0
        counter=counter+1;
        xk_store=xk_plus_sk_store;
        x_store(counter,:)=xk_store;
        f_of_xk_store=f_of_xk_plus_sk_store;
        if rho_k<eta_1
            delta=delta;
        elseif rho_k>=eta_1
            delta=min(gamma_1*delta,delta_max);
        end
        
        
        model_reg=QuadraticRegression([xk_store; gridpoints],EvalPoints);
        Constant=model_reg.Constant;
        Linear=model_reg.Linear;
        Squared=model_reg.Squared;
        polynomial_k=@(parameters) Constant+parameters*Linear+(1/2)*parameters*Squared*parameters';
        xk_plus_sk_store=fmincon(polynomial_k, xk_store, [], [] ,[], [], grid_lb, grid_ub,[], options);
    end
    x_store(counter,:)=xk_store
    tol=max(max(abs(xk_store-xk_plus_sk_store)),f_of_xk_plus_sk_store-f_of_xk_store);
    
    
end

out=x_store;
% matlabpool close
    
    
