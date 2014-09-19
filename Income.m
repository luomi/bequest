function y = Income(age, gender, category)

vec_constant    = [4.037 9.296 9.842 10.55 11.46]; 
vec_age2        = [-0.000677 -1.43e-05 -1.06e-05 4.21e-05 2.22e-05]; 
vec_age         = [0.114 0.00285 0.000911 -0.00636 -0.00616];
vec_gender      = [-0.0563 -0.108 -0.0524 -0.0261 0.0672] ;
vec_genderage   = [-1.55e-05 0.00145 0.000644 0.000430 -0.000273];

b_constant      = vec_constant(category);
b_age           = vec_age(category);
b_age2          = vec_age2(category);
b_gender        = vec_gender(category);
b_genderage     = vec_genderage(category);

l_inc           = b_constant+b_age*age+b_age2*age^2+b_gender*(gender==1)+b_genderage*(gender==1)*age;
y               = exp(l_inc)/1e3;