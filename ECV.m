function cont_val=ECV(a, state, a_grid,Value1,Value2,Value3,Value4, prob1_h, prob2_h,prob3_h,prob4_h, Health_prob, betta)


prob=Health_prob;

prob=prob(state,:);
[n,~]=size(Value1);

Value_int=zeros(n,4);
Value_int(:,1)=Value1*prob1_h; 
Value_int(:,2)=Value2*prob2_h;  
Value_int(:,3)=Value3*prob3_h; 
Value_int(:,4)=Value4*prob4_h; 

value_mat=[interp1(a_grid, Value_int(:,1), a,'cubic'),interp1(a_grid, Value_int(:,2), a,'cubic'),...
    interp1(a_grid, Value_int(:,3), a,'cubic'), interp1(a_grid, Value_int(:,4), a,'cubic')];

cont_val=betta*prob*value_mat';

    
