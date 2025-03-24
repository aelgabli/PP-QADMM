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
rho = 0.1; 
delta = 0.1;
epsilon = 0.1;
c1 = 15; % 1 or 2 or 5 (based on the ready-preprocessed dataset normalization)
bitsToSend = 3; % this is b, the # of bits to represent each model dimension (# of bits per sample)


num_iter = 40000;
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

acc = 1e-20;
transmissionTime = 1e-3; % this is tau = 1 ms which the upload/download transmission time

sigma = 50;%2*c1*sqrt(2*log(1.25/delta))/(noSamples*epsilon*rho);

% RUN THE CODE FOR
%----------------

% PS-ADMM
[obj_ADMM, loss_ADMM] = standard_ADMM ...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, acc);

% ADMM+Quantization (primal)
[obj_ADMM_w_Qnt, loss_ADMM_w_Qnt] = ADMM_w_Qnt...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend, acc);


 
    
% PP-QADMM
[obj_PPQADMM, loss_PPQADMM] = PPQADMM...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend, sigma, acc);

sigma = 100;
[obj_PPQADMM2, loss_PPQADMM2] = PPQADMM...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend, sigma, acc);

sigma = 200;
[obj_PPQADMM3, loss_PPQADMM3] = PPQADMM...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend, sigma, acc);

Iteration = 1:num_iter;
figure
semilogy(Iteration, loss_ADMM, Iteration, loss_ADMM_w_Qnt,...
    Iteration, loss_PPQADMM,Iteration, loss_PPQADMM2,Iteration, loss_PPQADMM3, 'LineWidth', 1.25)
grid on;
xlabel('Iteration')
ylabel('Loss')
legend('standard ADMM','QADMM'...
    ,'PP-QADMM, \sigma=5', 'PP-QADMM, \sigma=50','PP-QADMM, \sigma=100');


 %% Alert me after finishing the code
Data = load('splat.mat');  % handel   chirp   gong   train  splat
sound(Data.y, Data.Fs)

Total_time_in_minutes = toc/60

% End of code, thanks for watching :)

% save results_all_variousSigma.mat Iteration loss_ADMM loss_ADMM_w_Qnt...
%     loss_PPQADMM loss_PPQADMM2 loss_PPQADMM3...
%     no_workers noSamples

save results_100workers_variousSigma.mat Iteration loss_ADMM loss_ADMM_w_Qnt...
    loss_PPQADMM loss_PPQADMM2 loss_PPQADMM3...
    no_workers noSamples

% c1=15;
% delta=0.1;
% rho=0.1;
% sigma=200;
% epsilon=2*c1*sqrt(2*log(1.25/delta))/(noSamples*sigma*rho)
% 

