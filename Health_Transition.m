function prob=Health_Transition(a,gender)
a       = a-48; % transition starts at age 49 to 110
load('result_ADLH_all.mat','MK');
prob    = rot90([1,0,0,0;MK(:,:,gender,a)],2);

% prob    = [MK(:,:,1,a);0,0,0,1]; % male
% prob    = [MK(:,:,2,a);0,0,0,1]; % female

% P(1,1)=.963945;
% P(1,2)=.033547;
% P(1,3)=.000020;
% P(1,4)=1-sum(P(1,1:3));
% 
% P(2,1)=.336005;
% P(2,2)=.560655;
% P(2,3)=.065959;
% P(2,4)=1-sum(P(2,1:3));
% 
% P(3,1)=.024812;
% P(3,2)=.136231;
% P(3,3)=.746274;
% P(3,4)=1-sum(P(3,1:3));
% 
% 
% c1=.001441;
% c2=.8966;
% c3=.5643; 
% e=1.5;
% 
% AgeMat=eye(4,4);
% 
% AgeMat(1,1)=1-c1*a^e;
% AgeMat(1,2)=c1*a^e*(c2*c3)/(1+c2+c3*c2);
% AgeMat(1,3)=c1*a^e*(c3)/(1+c2+c3*c2);
% AgeMat(1,4)=c1*a^e*(1)/(1+c2+c3*c2);
% 
% 
% AgeMat(2,2)=1-c1*a^e;
% AgeMat(2,3)=c1*a^e*(c2)/(1+c2);
% AgeMat(2,4)=c1*a^e*(1)/(1+c2);
% 
% 
% AgeMat(3,3)=1-c1*a^e;
% AgeMat(3,4)=c1*a^e;
% 
% prob=P*AgeMat;