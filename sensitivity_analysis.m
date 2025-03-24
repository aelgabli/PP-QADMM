clc
% close all
clear all
tic % start timer


% READ THE DATA
%----------------
XX = readtable("data/x_train_MinMax_Normalized.csv"); % read the training data
XX = table2array(XX); % change the data type from table to double to allow manipulating the data
% csvread('data\x_train_MinMax_Normalized.csv',1,0);

YY = readtable("data/y_train_MinMax_Normalized.csv"); 
YY = table2array(YY); 

% Specify the normalized dataset type used to pass it to the graph title 
%-----------------------------------------------------------------------
 Dataset_type = 'MinMax-MaxNorm';
%Dataset_type = 'Zscore-MaxNorm';
% Dataset_type = 'MinMax';  
% Dataset_type = 'Zscore';

 

% DEFINE THE NEEDED PARAMETERS
%------------------------------
rho = 1; 
delta = 0.1;
epsilon = 0.1;
c1 = 15; % 1 or 2 or 5 (based on the ready-preprocessed dataset normalization)
bitsToSend = 3; % this is b, the # of bits to represent each model dimension (# of bits per sample)

%numWorkersArray=[10, 100, 500, 1000];

num_iter = 40000;

for jj=1:10

no_workers = 100;
%sigma =10;
num_feature = size(XX,2);
noSamples = floor(size(XX,1)/no_workers);
total_num_samples = noSamples * no_workers;

XX = XX(1:total_num_samples,:);
YY = YY(1:total_num_samples);

% solve for the optimal solution analytically via the least squares
% approach using ALL data (assuming centralized system that has all the data at its PS)
[w_optimal, obj0] = opt_sol_closedForm(XX,YY); 

acc = 1e-10;
transmissionTime = 1e-3; % this is tau = 1 ms which the upload/download transmission time

sigma = 2*c1*sqrt(2*log(1.25/delta))/(noSamples*epsilon*rho); 
%sigma = 2*c1*sqrt(2*log(1.25/delta))/(epsilon*rho); 

% RUN THE CODE FOR
%----------------

% PS-ADMM
[obj_ADMM, loss_ADMM, iter_admm] = standard_ADMM ...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, acc);

% ADMM+Quantization (primal)
[obj_ADMM_w_Qnt, loss_ADMM_w_Qnt, iter_qadmm] = ADMM_w_Qnt...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend, acc);



    
% PP-QADMM
[obj_PPQADMM, loss_PPQADMM, iter_ppqadmm] = PPQADMM...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend, sigma, acc);




d = 13; % the model dimension
b = 3; % number of bits per sample
w = zeros(1,d);
number_of_bits_to_Send = 64 + length(w)*b; 
%no_workers = 10;

total_bits_to_send_ours = zeros(1,iter_ppqadmm); % initialize
total_bits_to_send_QADMM = zeros(1,iter_qadmm); % initialize
total_bits_to_send_ADMM = zeros(1,iter_admm); % initialize



for i = 1:iter_ppqadmm

     if i == 1
              total_bits_to_send_ours(i) = number_of_bits_to_Send * (no_workers);
              
     else
         total_bits_to_send_ours(i) = total_bits_to_send_ours(i-1) + (number_of_bits_to_Send * (no_workers));
     

     end
end



for i = 1:iter_qadmm

     if i == 1
              total_bits_to_send_QADMM(i) = number_of_bits_to_Send * (no_workers);
     else
         total_bits_to_send_QADMM(i) = total_bits_to_send_QADMM(i-1) + (number_of_bits_to_Send * (no_workers));
        

     end
end



for i = 1:iter_admm

     if i == 1
              total_bits_to_send_ADMM(i) = length(w)*64* (no_workers);
     else
         total_bits_to_send_ADMM(i) = total_bits_to_send_ADMM(i-1) + length(w)*64* (no_workers);
         

     end
end




propQADMM(jj)=total_bits_to_send_QADMM(iter_qadmm)/total_bits_to_send_ADMM(iter_admm);
propPPQADMM(jj)=total_bits_to_send_ours(iter_ppqadmm)/total_bits_to_send_ADMM(iter_admm);
printoutMsg = sprintf('QADMM bit proportion is: %d .',propQADMM(jj));
disp(printoutMsg)
printoutMsg = sprintf('PPQADMM bit proportion is: %d .',propPPQADMM(jj));
disp(printoutMsg)





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

% X=100;
% Y=100;

X=100;
Y=100;

k = 3;

% Allocate the PS in the center of the squared region
x_PS = X/2;
y_PS = Y/2;

N0=1E-6; % The noise PSD N0 in Watt/Hz

% Initialize
d_square_central = zeros(1,num_workers);
P_central_uplink = zeros(1,num_workers);
E_total_qadmm = zeros(1,iter_qadmm);
E_total_ppqadmm = zeros(1,iter_ppqadmm);
E_total_admm = zeros(1,iter_admm);

% generate uniformly distributed (x,y) coordinates in the range feom 0 --> (X,Y)
for n=1:num_workers
  x(n)=X*rand;
  y(n)=Y*rand;
end


for n=1:num_workers
    d_square_central(n)=(x(n)-x_PS)^2+(y(n)-y_PS)^2; % find D^2 from the PS

    P_central_uplink_Q(n)= d_square_central(n)*N0*Band*(2^(2*Rate_Q/Band)-1); % P

    P_central_uplink_nonQ(n)= d_square_central(n)*N0*Band*(2^(2*Rate_nonQ/Band)-1); % P

end


for i = 1:iter_ppqadmm
        if i == 1
            E_total_ppqadmm(i) = sum(P_central_uplink_Q)*tau*i;
        
        else
            E_total_ppqadmm(i) = E_total_ppqadmm(i-1) + sum(P_central_uplink_Q)*tau*i;

        end
end



for i = 1:iter_qadmm
        if i == 1
            E_total_qadmm(i) = sum(P_central_uplink_Q)*tau*i;
        
        else
            E_total_qadmm(i) = E_total_qadmm(i-1) + sum(P_central_uplink_Q)*tau*i;

        end
end



for i = 1:iter_admm
        if i == 1
            E_total_admm(i) = sum(P_central_uplink_nonQ)*tau*i;
        
        else
            E_total_admm(i) = E_total_admm(i-1) + sum(P_central_uplink_nonQ)*tau*i;

        end
end




propQADMM_en(jj)=E_total_qadmm(iter_qadmm)/E_total_admm(iter_admm);
propPPQADMM_en(jj)=E_total_ppqadmm(iter_ppqadmm)/E_total_admm(iter_admm);
printoutMsg = sprintf('QADMM energy proportion is: %d .',propQADMM_en(jj));
disp(printoutMsg)
printoutMsg = sprintf('PPQADMM energy proportion is: %d .',propPPQADMM_en(jj));
disp(printoutMsg)



end


avg_propQADMM=sum(propQADMM)/10;
avg_propPPQADMM=sum(propPPQADMM)/10;
printoutMsg = sprintf('QADMM N-Bits: %d .',avg_propQADMM);
disp(printoutMsg)
printoutMsg = sprintf('PPQADMM N-Bits: %d .',avg_propPPQADMM);
disp(printoutMsg)


avg_propQADMM_en=sum(propQADMM_en)/10;
avg_propPPQADMM_en=sum(propPPQADMM_en)/10;
printoutMsg = sprintf('QADMM N-Energy: %d .',avg_propQADMM_en);
disp(printoutMsg)
printoutMsg = sprintf('PPQADMM N-Energy: %d .',avg_propPPQADMM_en);
disp(printoutMsg)

 %% Alert me after finishing the code
Data = load('splat.mat');  % handel   chirp   gong   train  splat
sound(Data.y, Data.Fs)

Total_time_in_minutes = toc/60

% End of code, thanks for watching :)

% save workers100.mat loss_ADMM loss_ADMM_w_Qnt ...
%     loss_PPQADMM no_workers noSamples

%save eff_noWorkers.mat propQADMM_en propPPQADMM_en propPPQADMM propQADMM



