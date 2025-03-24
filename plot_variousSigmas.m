%% PLOT Multiple figures

clear
close all

% load results_all_variousSigma.mat Iteration loss_ADMM loss_ADMM_w_Qnt...
%     loss_PPQADMM loss_PPQADMM2 loss_PPQADMM3...
%     no_workers noSamples

load results_100workers_variousSigma.mat Iteration loss_ADMM loss_ADMM_w_Qnt...
    loss_PPQADMM loss_PPQADMM2 loss_PPQADMM3...
    no_workers noSamples

Avg_loss_ADMM = loss_ADMM / (noSamples);
Avg_loss_ADMM_w_Qnt = loss_ADMM_w_Qnt / (noSamples);
Avg_loss_PPQADMM = loss_PPQADMM / (noSamples);
Avg_loss_PPQADMM2 = loss_PPQADMM2 / (noSamples);
Avg_loss_PPQADMM3 = loss_PPQADMM3 / (noSamples);




markerStep = 5000;


figure(1);
subplot(1,3,1)
semilogy(Avg_loss_ADMM,'MarkerIndices',1:markerStep:length(loss_ADMM),'LineWidth',2);
hold on
semilogy(Avg_loss_ADMM_w_Qnt,'--*','MarkerIndices',1:markerStep:length(loss_ADMM_w_Qnt),'LineWidth',2);
semilogy(Avg_loss_PPQADMM,':*k','MarkerIndices',1:markerStep:length(loss_PPQADMM),'LineWidth',2);
semilogy(Avg_loss_PPQADMM2,'--s','MarkerIndices',1:markerStep:length(loss_PPQADMM2),'LineWidth',2);
semilogy(Avg_loss_PPQADMM3,'--^','MarkerIndices',1:markerStep:length(loss_PPQADMM3),'LineWidth',2);
xlabel({'Number of Iterations';'(a)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Optimality Gap','fontsize',16,'fontname','Times New Roman')
legend('standard ADMM','QADMM'...
    ,'PP-QADMM, \sigma=50', 'PP-QADMM, \sigma=100','PP-QADMM, \sigma=200');
%title('\epsilon=0.1, \delta=0.1')
set(gca,'fontsize',14,'fontweight','bold');


d = 13; % the model dimension
b = 3; % number of bits per sample
w = zeros(1,d);
number_of_bits_to_Send = 64 + length(w)*b; 
no_workers = 10;
max = 40000;%1e6; % this is the max number of iteration to plot versus

total_bits_to_send_ours = zeros(1,max); % initialize
total_bits_to_send_QADMM = zeros(1,max); % initialize
total_bits_to_send_ADMM = zeros(1,max); % initialize


for i = 1:max


% Calculate the total cumulative number of Tx'd bits till itaeration i
     %---------------------------------------------------------------------
     if i == 1
              total_bits_to_send_ours(i) = number_of_bits_to_Send * (no_workers);
              total_bits_to_send_QADMM(i) = number_of_bits_to_Send * (no_workers);
              total_bits_to_send_ADMM(i) = length(w)*64* (no_workers);
     else
         total_bits_to_send_ours(i) = total_bits_to_send_ours(i-1) + (number_of_bits_to_Send * (no_workers));
         total_bits_to_send_QADMM(i) = total_bits_to_send_QADMM(i-1) + (number_of_bits_to_Send * (no_workers));
         total_bits_to_send_ADMM(i) = total_bits_to_send_ADMM(i-1) + length(w)*64* (no_workers);

     end

end



figure(1);
subplot(1,3,2)

markerStep = 10000;

semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM,'MarkerIndices',1:markerStep:length(loss_ADMM),'LineWidth',2);
hold on
semilogy(total_bits_to_send_QADMM,Avg_loss_ADMM_w_Qnt,'--*','MarkerIndices',1:markerStep:length(loss_ADMM_w_Qnt),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP,'--s','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP2,'--^','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP2),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP3,':D','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP3),'LineWidth',2);
semilogy(total_bits_to_send_ours,Avg_loss_PPQADMM,':*k','MarkerIndices',1:markerStep:length(loss_PPQADMM),'LineWidth',2);
semilogy(total_bits_to_send_ours,Avg_loss_PPQADMM2,'--s','MarkerIndices',1:markerStep:length(Avg_loss_PPQADMM2),'LineWidth',2);
semilogy(total_bits_to_send_ours,Avg_loss_PPQADMM3,'--^','MarkerIndices',1:markerStep:length(Avg_loss_PPQADMM3),'LineWidth',2);
xlabel({'Total Number of Transmitted Bits';'(b)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Optimality Gap','fontsize',16,'fontname','Times New Roman')
% legend('standard ADMM','QADMM','DPADMM'...
%     ,'DPADMM, \rho=5','DPADMM, \rho=50','PPQADMM');
xlim([0 2E8])

legend('standard ADMM','QADMM'...
    ,'PP-QADMM, \sigma=50', 'PP-QADMM, \sigma=100','PP-QADMM, \sigma=200');
%title('\epsilon=0.1, \delta=0.1')
set(gca,'fontsize',14,'fontweight','bold');


% This script to plot the loss vs the total consumed energy 

d = 13; % the model dimension
b = 3;
w = zeros(1,d);
num_workers = 10;
BW = 1e6; % 1 MHz
Band = BW / num_workers;
tau = 1e-7; % 0.1 msec
N0 = 1e-6; % The noise PSD N0 in Watt/Hz
Rate_Q = 64 + length(w)*b; % the rate for OURs and ADMM+Quant
Rate_nonQ = length(w)*64; % the rate for the remaining schemes w/o Quant
X=100;
Y=100;
max = 40000;
k = 3;

% Allocate the PS in the center of the squared region
x_PS = X/2;
y_PS = Y/2;

N0=1E-6; % The noise PSD N0 in Watt/Hz

% Initialize
d_square_central = zeros(1,num_workers);
P_central_uplink = zeros(1,num_workers);
E_total_Q = zeros(1,max);
E_total_nonQ = zeros(1,max);

% generate uniformly distributed (x,y) coordinates in the range feom 0 --> (X,Y)
for n=1:num_workers
  x(n)=X*rand;
  y(n)=Y*rand;
end


for n=1:num_workers
    d_square_central(n)=(x(n)-x_PS)^2+(y(n)-y_PS)^2; % find D^2 from the PS

    P_central_uplink_Q(n)= d_square_central(n)*N0*Band*(2^(2*Rate_Q/Band)-1); % P

    P_central_uplink_nonQ(n)= d_square_central(n)*N0*Band*(2^(2*Rate_nonQ/Band)-1); % P

    % P_central_downlink(n)= num_workers/2*Band*d_square_central(n)*N0*(2^(2*Rate/(num_workers*Band))-1);
end

for i = 1:max
        if i == 1
            E_total_Q(i) = sum(P_central_uplink_Q)*tau*i;

            E_total_nonQ(i) = sum(P_central_uplink_nonQ)*tau*i;
        
        else
            E_total_Q(i) = E_total_Q(i-1) + sum(P_central_uplink_Q)*tau*i;

            E_total_nonQ(i) =  E_total_nonQ(i-1)  + sum(P_central_uplink_nonQ)*tau*i;
        end
end


figure(1);
subplot(1,3,3)

markerStep = 10000;

semilogy(E_total_nonQ,Avg_loss_ADMM,'MarkerIndices',1:markerStep:length(loss_ADMM),'LineWidth',2);
hold on
semilogy(E_total_Q,Avg_loss_ADMM_w_Qnt,'--*','MarkerIndices',1:markerStep:length(loss_ADMM_w_Qnt),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP,'--s','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP2,'--^','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP2),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP3,':D','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP3),'LineWidth',2);
semilogy(E_total_Q,Avg_loss_PPQADMM,':*k','MarkerIndices',1:markerStep:length(loss_PPQADMM),'LineWidth',2);
semilogy(E_total_Q,Avg_loss_PPQADMM2,'--s','MarkerIndices',1:markerStep:length(Avg_loss_PPQADMM2),'LineWidth',2);
semilogy(E_total_Q,Avg_loss_PPQADMM3,'--^','MarkerIndices',1:markerStep:length(Avg_loss_PPQADMM3),'LineWidth',2);
xlabel({'Total Energy';'(c)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Optimality Gap','fontsize',16,'fontname','Times New Roman')
xlim([0 600])
legend('standard ADMM','QADMM'...
    ,'PP-QADMM, \sigma=50', 'PP-QADMM, \sigma=100','PP-QADMM, \sigma=200');
%title('\epsilon=0.1, \delta=0.1')
set(gca,'fontsize',14,'fontweight','bold');


