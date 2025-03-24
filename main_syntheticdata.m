clc
% close all
clear all
tic % start timer


% READ THE DATA
%----------------

%XX_init = readtable("data/zscore_normalized_with_bias_synthetic_data.csv");
XX_init = readtable("data/fully_normalized_with_bias_feature.csv"); % read the training data

XX = table2array(XX_init(1:end,1:50)); % change the data type from table to double to allow manipulating the data
YY = table2array(XX_init(1:end, 51));

YY=YY+randn(size(YY,1),1);


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
c1 = 52; % 1 or 2 or 5 (based on the ready-preprocessed dataset normalization)
bitsToSend = 3; % this is b, the # of bits to represent each model dimension (# of bits per sample)



num_iter = 3000;
no_workers = 100;
%sigma =10;
num_feature = size(XX,2);
noSamples = floor(size(XX,1)/no_workers);
total_num_samples = noSamples * no_workers;

XX = XX(1:total_num_samples,:);
YY = YY(1:total_num_samples);

% solve for the optimal solution analytically via the least squares
[w_optimal, obj0] = opt_sol_closedForm(XX,YY); 

acc = 1e-30;
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






%% PLOT
Iteration = 1:num_iter;


figure
semilogy(1:iter_admm, loss_ADMM,...
     1:iter_qadmm, loss_ADMM_w_Qnt,...
    1:iter_ppqadmm, loss_PPQADMM,...
    'LineWidth', 1.25)
grid on;
xlabel('Iteration')
ylabel('Loss')
legend ('ADMM','QADMM','PPQADMM')


 %% Alert me after finishing the code
Data = load('splat.mat');  % handel   chirp   gong   train  splat
sound(Data.y, Data.Fs)

Total_time_in_minutes = toc/60

% End of code, thanks for watching :)

% save results_100workers_synthetic.mat Iteration loss_ADMM loss_ADMM_w_Qnt...
%     loss_PPQADMM no_workers noSamples

