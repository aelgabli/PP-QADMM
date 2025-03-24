clc
clear all
tic % start timer


% READ THE DATA
%----------------
XX = readtable("data/x_train_MinMax_Normalized.csv"); % read the training data
XX = table2array(XX); % change the data type from table to double to allow manipulating the data


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
[w_optimal, obj0] = opt_sol_closedForm(XX,YY); 

acc = 1e-20;
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


% ADMM+DP
[obj_ADMM_w_DP, loss_ADMM_w_DP] = ADMM_w_DP...
           (XX,YY, rho, delta, epsilon, no_workers, num_feature, noSamples, num_iter, obj0,c1);



c1 = 5;
 
XX = readtable("data/X_train_MinMax_MaxNormalized (c1=5).csv"); 
XX = table2array(XX); % change the data type from table to double to allow manipulating the data

YY = readtable("data/y_train_MinMax_Normalized.csv"); 
YY = table2array(YY); 

XX = XX(1:total_num_samples,:);
YY = YY(1:total_num_samples);

% solve for the optimal solution analytically via the least squares
[w_optimal, obj0] = opt_sol_closedForm(XX,YY); 


rho=5;
% 
% ADMM+DP
[obj_ADMM_w_DP2, loss_ADMM_w_DP2] = ADMM_w_DP...
           (XX,YY, rho, delta, epsilon, no_workers, num_feature, noSamples, num_iter, obj0,c1);

rho=50;
% 
% ADMM+DP
[obj_ADMM_w_DP3, loss_ADMM_w_DP3] = ADMM_w_DP...
           (XX,YY, rho, delta, epsilon, no_workers, num_feature, noSamples, num_iter, obj0,c1);



%% PLOT
Iteration = 1:num_iter;

figure
semilogy(Iteration, loss_ADMM, Iteration, loss_ADMM_w_Qnt, Iteration,...
    loss_ADMM_w_DP, Iteration, loss_ADMM_w_DP2,...
    Iteration, loss_ADMM_w_DP3,...
    Iteration, loss_PPQADMM, 'LineWidth', 1.25)
grid on;
xlabel('Iteration')
ylabel('Loss')
legend ('ADMM','QADMM', 'DP-ADMM','DP-ADMM,\rho=5', 'DP-ADMM,\rho=50','PPQADMM')

% title('PS-ADMM Loss vs. Iterations (W/o the intercept)')

% title(['[# of Workers= ',num2str(no_workers),', \epsilon=',num2str(epsilon), ...
%     ', \delta=', num2str(delta),', \rho=',num2str(rho),', Qnt-bits= ',num2str(bitsToSend),...
%     ', c1=', num2str(c1),', Data: ',Dataset_type,']'], 'FontSize',15)


 %% Alert me after finishing the code
Data = load('splat.mat');  % handel   chirp   gong   train  splat
sound(Data.y, Data.Fs)

Total_time_in_minutes = toc/60

% End of code, thanks for watching :)

% save results_100workers.mat Iteration loss_ADMM loss_ADMM_w_Qnt loss_ADMM_w_DP...
%     loss_PPQADMM loss_ADMM_w_DP2 loss_ADMM_w_DP3...
%     no_workers noSamples

