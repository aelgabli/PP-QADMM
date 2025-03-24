%% PLOT Multiple figures

clear
close all

% load results_all.mat Iteration loss_ADMM loss_ADMM_w_Qnt loss_ADMM_w_DP...
%     loss_PPQADMM loss_ADMM_w_DP2 loss_ADMM_w_DP3...
%     no_workers noSamples

load results_100workers_noniid.mat Iteration loss_ADMM loss_ADMM_w_Qnt loss_ADMM_w_DP...
    loss_PPQADMM loss_ADMM_w_DP2 loss_ADMM_w_DP3...
    no_workers noSamples


Avg_loss_ADMM = loss_ADMM ;
Avg_loss_ADMM_w_Qnt = loss_ADMM_w_Qnt;
Avg_loss_ADMM_w_DP = loss_ADMM_w_DP ;
Avg_loss_ADMM_w_DP2 = loss_ADMM_w_DP2;
Avg_loss_ADMM_w_DP3 = loss_ADMM_w_DP3;
Avg_loss_PPQADMM = loss_PPQADMM;




markerStep = 5000;

Avg_loss_HE = Avg_loss_ADMM;
Avg_loss_MPC = Avg_loss_ADMM;

figure(1);
subplot(1,3,1)
semilogy(Avg_loss_ADMM,'MarkerIndices',1:markerStep:length(loss_ADMM),'LineWidth',2);
hold on
semilogy(Avg_loss_HE,':','MarkerIndices',1:markerStep:length(Avg_loss_HE),'LineWidth',2);
semilogy(Avg_loss_MPC,'--','MarkerIndices',1:markerStep:length(Avg_loss_MPC),'LineWidth',2);
semilogy(Avg_loss_ADMM_w_Qnt,'--*','MarkerIndices',1:markerStep:length(loss_ADMM_w_Qnt),'LineWidth',2);
semilogy(Avg_loss_ADMM_w_DP,'--s','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP),'LineWidth',2);
semilogy(Avg_loss_ADMM_w_DP2,'--^','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP2),'LineWidth',2);
semilogy(Avg_loss_ADMM_w_DP3,':D','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP3),'LineWidth',2);
semilogy(Avg_loss_PPQADMM,':*k','MarkerIndices',1:markerStep:length(loss_PPQADMM),'LineWidth',2);
xlabel({'Number of Iterations';'(a)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Optimality Gap','fontsize',16,'fontname','Times New Roman')
legend('standard ADMM, \rho=0.1','HE-ADMM, \rho=0.1','MPC-ADMM, \rho=0.1','QADMM, \rho=0.1','DPADMM, \rho=0.1'...
    ,'DPADMM, \rho=5','DPADMM, \rho=50','PP-QADMM, \rho=0.1');
%title('\epsilon=0.1, \delta=0.1')
set(gca,'fontsize',14,'fontweight','bold');


d = 13; % the model dimension
b = 3; % number of bits per sample
w = zeros(1,d);
number_of_bits_to_Send = 64 + length(w)*b; 
%no_workers = 10;
max = 40000;%1e6; % this is the max number of iteration to plot versus

total_bits_to_send_ours = zeros(1,max); % initialize
total_bits_to_send_QADMM = zeros(1,max); % initialize
total_bits_to_send_ADMM = zeros(1,max); % initialize
total_bits_HE = zeros(1,max);
total_bits_MPC = zeros(1,max);


flagHE=0;
flagQPPADMM=0;
flagMPC=0;
flagADMM = 0;
for i = 1:max

     if i == 1
              total_bits_to_send_ours(i) = number_of_bits_to_Send * (no_workers);
              total_bits_to_send_QADMM(i) = number_of_bits_to_Send * (no_workers);
              total_bits_to_send_ADMM(i) = length(w)*64* (no_workers);
              total_bits_HE(i) = d*4096*no_workers;
              total_bits_MPC(i) = d*64*no_workers*no_workers;
     else
         total_bits_to_send_ours(i) = total_bits_to_send_ours(i-1) + (number_of_bits_to_Send * (no_workers));
         total_bits_to_send_QADMM(i) = total_bits_to_send_QADMM(i-1) + (number_of_bits_to_Send * (no_workers));
         total_bits_to_send_ADMM(i) = total_bits_to_send_ADMM(i-1) + length(w)*64* (no_workers);
         total_bits_HE(i) = total_bits_HE(i-1)+d*4096*no_workers;
         total_bits_MPC(i) = total_bits_MPC(i-1)+d*64*no_workers*no_workers;

     end

     if (Avg_loss_ADMM(i) <= 1E-10 && flagADMM==0 )
         bitsADMM=total_bits_to_send_ADMM(i);
        XX = sprintf('HE bits are: %d .',total_bits_to_send_ADMM(i));
        disp(XX)
        flagADMM =1;
     end

    if (Avg_loss_HE(i) <= 1E-10 && flagHE==0 )
        bitsHE=total_bits_HE(i);
        XX = sprintf('HE bits are: %d .',total_bits_HE(i));
        disp(XX)
        flagHE =1;
    end

    if (Avg_loss_MPC(i) <= 1E-10 && flagMPC==0 )
        bitsMPC=total_bits_MPC(i);
        XX = sprintf('MPC bits are: %d .',total_bits_MPC(i));
        disp(XX)
        flagMPC =1;
    end

    if (Avg_loss_PPQADMM(i) <= 1E-10 && flagQPPADMM==0 )
        bitsPPQADMM=total_bits_to_send_ours(i);
        XX = sprintf('QPPADMM bits are: %d .',total_bits_to_send_ours(i));
        disp(XX)
        flagQPPADMM =1;
    end


end

propHE=bitsHE/bitsADMM;
propMPC=bitsMPC/bitsADMM;
propPPQADMM=bitsPPQADMM/bitsADMM;
XX = sprintf('HE proportion is: %d .',propHE);
disp(XX)
XX = sprintf('MPC proportion is: %d .',propMPC);
disp(XX)
XX = sprintf('QPPADMM proportion is: %d .',propPPQADMM);
disp(XX)


figure(1);
subplot(1,3,2)

markerStep = 10000;
semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM,'MarkerIndices',1:markerStep:length(loss_ADMM),'LineWidth',2);
hold on
semilogy(total_bits_HE,Avg_loss_HE,':','MarkerIndices',1:markerStep:length(Avg_loss_HE),'LineWidth',2);
semilogy(total_bits_MPC,Avg_loss_MPC,'--','MarkerIndices',1:markerStep:length(Avg_loss_MPC),'LineWidth',2);
semilogy(total_bits_to_send_QADMM,Avg_loss_ADMM_w_Qnt,'--*','MarkerIndices',1:markerStep:length(loss_ADMM_w_Qnt),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP,'--s','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP2,'--^','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP2),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP3,':D','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP3),'LineWidth',2);
semilogy(total_bits_to_send_ours,Avg_loss_PPQADMM,':*k','MarkerIndices',1:markerStep:length(loss_PPQADMM),'LineWidth',2);
xlabel({'Total Number of Transmitted Bits';'(b)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Optimality Gap','fontsize',16,'fontname','Times New Roman')
% legend('standard ADMM','QADMM','DPADMM'...
%     ,'DPADMM, \rho=5','DPADMM, \rho=50','PPQADMM');
xlim([0 2E9])
legend('standard ADMM, \rho=0.1','HE-ADMM, \rho=0.1','MPC-ADMM, \rho=0.1','QADMM, \rho=0.1','PP-QADMM, \rho=0.1');
%title('\epsilon=0.1, \delta=0.1')
set(gca,'fontsize',14,'fontweight','bold');


% This script to plot the loss vs the total consumed energy 

d = 13; % the model dimension
b = 3;
w = zeros(1,d);
num_workers = no_workers;
BW = 1e6; % 1 MHz
Band = BW / num_workers;
tau = 1e-7; % 0.1 msec
N0 = 1e-6; % The noise PSD N0 in Watt/Hz
Rate_Q = 64 + length(w)*b; % the rate for OURs and ADMM+Quant
Rate_nonQ = length(w)*64; % the rate for the remaining schemes w/o Quant
Rate_HE=length(w)*4096;
Rate_MPC=length(w)*64*num_workers;
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

    P_central_uplink_HE(n)= d_square_central(n)*N0*Band*(2^(2*Rate_HE/Band)-1);
    P_central_uplink_MPC(n)= d_square_central(n)*N0*Band*(2^(2*Rate_MPC/Band)-1);

    % P_central_downlink(n)= num_workers/2*Band*d_square_central(n)*N0*(2^(2*Rate/(num_workers*Band))-1);
end

flagHE=0;
flagQPPADMM=0;
flagMPC=0;
flagADMM = 0;
for i = 1:max
        if i == 1
            E_total_Q(i) = sum(P_central_uplink_Q)*tau*i;

            E_total_nonQ(i) = sum(P_central_uplink_nonQ)*tau*i;

            E_total_HE(i) = sum(P_central_uplink_HE)*tau*i;
            E_total_MPC(i) = sum(P_central_uplink_MPC)*tau*i;
        
        else
            E_total_Q(i) = E_total_Q(i-1) + sum(P_central_uplink_Q)*tau*i;

            E_total_nonQ(i) =  E_total_nonQ(i-1)  + sum(P_central_uplink_nonQ)*tau*i;
            E_total_HE(i) =  E_total_HE(i-1)  + sum(P_central_uplink_HE)*tau*i;
            E_total_MPC(i) =  E_total_MPC(i-1)  + sum(P_central_uplink_MPC)*tau*i;
        end

    if (Avg_loss_ADMM(i) <= 1E-10 && flagADMM==0 )
         energyADMM=E_total_nonQ(i);
        XX = sprintf('HE bits are: %d .',E_total_nonQ(i));
        disp(XX)
        flagADMM =1;
     end

    if (Avg_loss_HE(i) <= 1E-10 && flagHE==0 )
        energyHE = E_total_HE(i);
        XX = sprintf('HE energy is: %d .',E_total_HE(i));
        disp(XX)
        flagHE =1;
    end

    if (Avg_loss_MPC(i) <= 1E-10 && flagMPC==0 )
        energyMPC = E_total_MPC(i);
        XX = sprintf('MPC energy is: %d .',E_total_MPC(i));
        disp(XX)
        flagMPC =1;
    end

    if (Avg_loss_PPQADMM(i) <= 1E-10 && flagQPPADMM==0 )
        energyPPQADMM = E_total_Q(i);
        XX = sprintf('QPPADMM energy is: %d .',E_total_Q(i));
        disp(XX)
        flagQPPADMM =1;
    end

end

propHE=energyHE/energyADMM;
propMPC=energyMPC/energyADMM;
propPPQADMM=energyPPQADMM/energyADMM;
XX = sprintf('HE proportion is: %d .',propHE);
disp(XX)
XX = sprintf('MPC proportion is: %d .',propMPC);
disp(XX)
XX = sprintf('QPPADMM proportion is: %d .',propPPQADMM);
disp(XX)


figure(1);
subplot(1,3,3)

markerStep = 10000;

semilogy(E_total_nonQ,Avg_loss_ADMM,'MarkerIndices',1:markerStep:length(loss_ADMM),'LineWidth',2);
hold on
semilogy(E_total_HE,Avg_loss_HE,':','MarkerIndices',1:markerStep:length(Avg_loss_HE),'LineWidth',2);
semilogy(E_total_MPC,Avg_loss_MPC,'--','MarkerIndices',1:markerStep:length(Avg_loss_MPC),'LineWidth',2);
semilogy(E_total_Q,Avg_loss_ADMM_w_Qnt,'--*','MarkerIndices',1:markerStep:length(loss_ADMM_w_Qnt),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP,'--s','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP2,'--^','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP2),'LineWidth',2);
%semilogy(total_bits_to_send_ADMM,Avg_loss_ADMM_w_DP3,':D','MarkerIndices',1:markerStep:length(loss_ADMM_w_DP3),'LineWidth',2);
semilogy(E_total_Q,Avg_loss_PPQADMM,':*k','MarkerIndices',1:markerStep:length(loss_PPQADMM),'LineWidth',2);
xlabel({'Total Energy';'(c)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Optimality Gap','fontsize',16,'fontname','Times New Roman')
xlim([0 6000])
% legend('standard ADMM','QADMM','DPADMM'...
%     ,'DPADMM, \rho=5','DPADMM, \rho=50','PPQADMM');
legend('standard ADMM, \rho=0.1','HE-ADMM, \rho=0.1','MPC-ADMM, \rho=0.1','QADMM, \rho=0.1','PP-QADMM, \rho=0.1');
%title('\epsilon=0.1, \delta=0.1')
set(gca,'fontsize',14,'fontweight','bold');

