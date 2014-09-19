clc
clear

%% Data

data_trans    = csvread('core_trans_in.csv');
age     = data_trans(:,1);
trans   = data_trans(:,2)/1000;

age_mat = reshape(55:90,[3,12])';

age_mat_before=[zeros(1,3);age_mat(1:11,:)];
age_mat_next=[age_mat(2:12,:);zeros(1,3)];
age_mat=[age_mat_before,age_mat,age_mat_next];



trans_mom_25=zeros(12,1);
trans_mom_50=zeros(12,1);
trans_mom_75=zeros(12,1);

for ii=1:12
    trans_temp=trans(ismember(age, age_mat(ii,:)));
    cell_count=sum(ismember(age, age_mat(ii,:)));

 
%     trans_mom_10(ii)=prctile(trans_temp,10);
    trans_mom_25(ii)=prctile(trans_temp,25);
    trans_mom_50(ii)=prctile(trans_temp,50);
    trans_mom_75(ii)=prctile(trans_temp,75);
%     trans_mom_90(ii)=prctile(trans_temp,90);
end
Data_trans_mom=[trans_mom_25; trans_mom_50; trans_mom_75];

%% Plots
agevec  = linspace(56,89,12)';
plot(agevec,assets_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(25:36),'b+:','LineWidth',2);hold on;
legend('Model','Data','Location','Northwest');
plot(agevec,assets_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(13:24),'b+:','LineWidth',2);hold on;
plot(agevec,assets_moments(1:12),'ro-','LineWidth',2);hold on;
plot(agevec,Data_mom(1:12),'b+:','LineWidth',2);xlabel('Age');ylabel('Assets in thousands');hold off;

agevec  = linspace(56,89,12)';
plot(agevec,assets_moments(25:36),'ro-','LineWidth',2);hold on;
% plot(agevec,Data_mom(25:36),'b+:','LineWidth',2);hold on;
legend('Model','Data','Location','Northwest');
plot(agevec,assets_moments(13:24),'ro-','LineWidth',2);hold on;
% plot(agevec,Data_mom(13:24),'b+:','LineWidth',2);hold on;
plot(agevec,assets_moments(1:12),'ro-','LineWidth',2);hold on;
%plot(agevec,Data_mom(1:12),'b+:','LineWidth',2);
xlabel('Age');ylabel('Assets in thousands');hold off;

agevec  = linspace(56,89,12)';
% plot(agevec,cf_moments(25:36),'ro-','LineWidth',2);hold on;
plot(agevec,Data_trans_mom(25:36),'b+:','LineWidth',2);hold on;
% legend('Model','Data','Location','Northwest');
% plot(agevec,cf_moments(13:24),'ro-','LineWidth',2);hold on;
plot(agevec,Data_trans_mom(13:24),'b+:','LineWidth',2);hold on;
% plot(agevec,cf_moments(1:12),'ro-','LineWidth',2);hold on;
plot(agevec,Data_trans_mom(1:12),'b+:','LineWidth',2);
xlabel('Age');ylabel('Transfer in thousands');
hold off;

