function [obj_PPQADMM, loss_PPQADMM, Iter] = PPQADMM...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend, sigma, acc)
      

Iter= num_iter;

%lambda = zeros(num_feature,no_workers);
lambda = sigma*randn(num_feature,no_workers);

prev_out=zeros(num_feature,no_workers);
local_model=zeros(num_feature,no_workers);
quantized=zeros(num_feature,no_workers); % quantized model

out=zeros(num_feature,no_workers); % the local models of the workers
out_central=zeros(num_feature,1); % the global model of the parameter server

max_iter = num_iter;


     clc
     msg = ['Running Simulation for "ADMM+Quantization", hold on ... (2 out of 5)'];
     disp(msg)

 for i = 1:max_iter
     
     clc
     msg = ['Itaeration # ', num2str(i), ' for "ADMM+Quantization", hold on ... (2 out of 5) '];
     disp(msg)
     
     
     % to update the local model parameters (small thetas)
     for ii =1:no_workers

         term_1=rho*out_central;

         B1 = lambda(:,ii); 
         first = (ii-1)*noSamples+1; % specify the index for the first data sample for each user
         last = first+noSamples-1; % and also the last one.
        
        X=XX(first:last,1:num_feature); % Form the dataset for each worker (rows for different data samples, while columns for different features for each sample)
        Y=YY(first:last); % Form the corresponding target variable set

        x=((1/noSamples)*(X'*X)+rho*eye(num_feature))\((1/noSamples)*(X'*Y)-B1+term_1);
        % x=(X'*X+rho*eye(num_feature))\(X'*Y-B1+term_1); % x is the local model
        local_model(:,ii)=x;
        out(:,ii) =x+1/rho*B1;

         % quantized(:,ii) = x;
        % For mapping between variables below and in the fn itself : out == current
       [quantized(:,ii),number_of_bits_toSend]=stochasticQuantizer ...
                                        (quantized(:,ii),out(:,ii),prev_out(:,ii),bitsToSend);
        
     end
            prev_out = quantized;  % this is to update the previously quantized model

    % R
    % Delta
    
    out_central=zeros(num_feature,1);  % reset the global model to update it from the new local models

    % update the global model parameters (capital Theta)
    for ii =1:no_workers
       
        out_central = out_central +  1/(no_workers)*quantized(:,ii);
    end
     
         
    % update the dual variables as in Eq. (12)
    for ii=1:no_workers

      lambda(:,ii) = lambda(:,ii) + rho*(local_model(:,ii)-out_central);

    end
    

    
         
        final_obj = 0.5*norm(XX * out_central - YY)^2; % Note: this is the sum, we will compute the average in the plot file
        
        obj_PPQADMM(i)=final_obj;
        loss_PPQADMM(i)=abs(final_obj-obj0) ;    % obj0 is the optimal global model obtained by running the fn opt_sol_closedForm       
        
     if (loss_PPQADMM(i) <= acc)
         Iter = i;
        printoutMsg = sprintf('number of iterations PPQADMM: %d .',i);
        disp(printoutMsg)
            break;
        end
        
       

 end % end for i =1: maxiter   
    
     clc
     msg = ['Simulation finished for "ADMM+Quantization", moving forward ... '];
     disp(msg)

end % end function
     




