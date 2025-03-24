function [obj_ADMM_w_Qnt, loss_ADMM_w_Qnt] = ADMM_w_Qnt_noniid...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, bitsToSend)
      

Iter= num_iter;              
lambda = zeros(num_feature,no_workers);

prev_out=zeros(num_feature,no_workers); 
quantized=zeros(num_feature,no_workers); % quantized model

out=zeros(num_feature,no_workers); % the local models of the workers
out_central=zeros(num_feature,1); % the global model of the parameter server

max_iter = num_iter;


% com_cost(1)= 0;
% upload_slots(1)= 0;
% b=32;
% duration=1E-3; % this is tau = 1 ms which the upload/download transmission time

     clc
     msg = ['Running Simulation for "ADMM+Quantization", hold on ... (2 out of 5)'];
     disp(msg)

 for i = 1:max_iter
     % if(i > 1)
     %     com_cost(i)=com_cost(i-1);
     %     upload_slots(i)= upload_slots(i-1);
     % end
     clc
     msg = ['Itaeration # ', num2str(i), ' for "ADMM+Quantization", hold on ... (2 out of 5) '];
     disp(msg)
     
     
     % clc
     % msg = ['Current iteration: ', num2str(i)];
     % disp(msg)
     
     
     % to update the local model parameters (small thetas)
     for ii =1:no_workers

         term_1=rho*out_central;

         B1 = lambda(:,ii); % initialy this will be all-zeros, then it will be updated below.
         
         first = sum(noSamples(1:ii-1))+1; % specify the index for the first data sample for each user
         last = first+noSamples(ii)-1; % and also the last one.
        
        X=XX(first:last,1:num_feature); % Form the dataset for each worker (rows for different data samples, while columns for different features for each sample)
        Y=YY(first:last); % Form the corresponding target variable set

        x=((1/noSamples(ii))*(X'*X)+rho*eye(num_feature))\((1/noSamples(ii))*(X'*Y)-B1+term_1);
        % x=(X'*X+rho*eye(num_feature))\(X'*Y-B1+term_1); % x is the local model
        out(:,ii) =x;

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
       
        out_central = out_central +  1/(no_workers)*(quantized(:,ii)+lambda(:,ii)/rho);
    end
     
         
    % update the dual variables as in Eq. (12)
    for ii=1:no_workers

      lambda(:,ii) = lambda(:,ii) + rho*(quantized(:,ii)-out_central);

    end
    

    
         
        final_obj=0;
        for j=1:no_workers
            first = sum(noSamples(1:j-1))+1;
            last = first+noSamples(j)-1;
                
            X=XX(first:last,1:num_feature);
            Y=YY(first:last);
            final_obj =final_obj+1/noSamples(j)*0.5*norm(X * out_central - Y)^2;
        end
        
        obj_ADMM_w_Qnt(i)=final_obj;
        loss_ADMM_w_Qnt(i)=abs(final_obj-obj0) ;    % obj0 is the optimal global model obtained by running the fn opt_sol_closedForm       
        
     
        
        % min_loss = loss_QADMM(i)

        
        
       

 end % end for i =1: maxiter   
    
     clc
     msg = ['Simulation finished for "ADMM+Quantization", moving forward ... '];
     disp(msg)

end % end function
     




