function val=g(parameter,e_states, a_grid, T,r,data , Data_mom, W)


f=@(q) Sim_Moments(q,e_states,a_grid, T,r, data);
% format long
disp(parameter)
% id=labindex;
% filename_param=['P:\nyu\working\codes\Structural_Model_Est\Save_Params\pval' num2str(id) '.txt'];
% filename_fval=['P:\nyu\working\codes\Structural_Model_Est\Save_Params\fval' num2str(id) '.txt'];

eval_mom=f(parameter);
val=(eval_mom-Data_mom)'*W*(eval_mom-Data_mom);

% dlmwrite(filename_param,parameter,  '-append','roffset',1,'delimiter',';','newline', 'pc');
% dlmwrite(filename_fval, val , '-append','roffset',1,'delimiter',';','newline', 'pc');