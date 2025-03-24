
function [obj_DPQADMM, loss_DPQADMM] = ADMM_w_DP...
    (XX,YY, rho, delta, epsilon, no_workers, num_feature, noSamples, num_iter, obj0,c1)
      


Iter= num_iter;              
lambda = zeros(num_feature,no_workers);

out=zeros(num_feature,no_workers); %  the local models of the workers
out_central=zeros(num_feature,1); %  the global model of the parameter server

max_iter = num_iter;



     clc
     msg = ['Running Simulation for "ADMM+DP", hold on ... (3 out of 5) '];
     disp(msg)


 for i = 1:max_iter
    
     
     clc
     msg = ['Itaeration # ', num2str(i), ' for "ADMM+DP", hold on ... (3 out of 5) '];
     disp(msg)
         
     
     % to update the local model parameters (small thetas)
     for ii =1:no_workers

         term_1=rho*out_central;

         B1 = lambda(:,ii); % initialy this will be all-zeros, then it will be updated below.
         
         first = (ii-1)*noSamples+1; % specify the index for the first data sample for each user
         last = first+noSamples-1; % and also the last one.
        
         % Form the dataset for each worker (rows for different data samples, while columns for different features for each sample)
        X=XX(first:last,1:num_feature); 
        Y=YY(first:last); % Form the corresponding target variable set

        % x=(X'*X+rho*eye(num_feature,num_feature))\(X'*Y-B1+term_1);
        x=((1/noSamples)*(X'*X)+rho*eye(num_feature))\((1/noSamples)*(X'*Y)-B1+term_1);

        out(:,ii) =x;

        %-----------------------------------------
        % Adding Gaussian noise to the local model 
        %-----------------------------------------

        % define the needed parameters
        % define the std. dev. to control the noise variance
        sigma = 2*c1*sqrt(2*log(1.25/delta))/(noSamples*epsilon*rho); 
        % sigma = 2*c1*sqrt(2*log(1.25/delta))/(epsilon*rho);  

        Gauss_noise = sigma*randn(num_feature,1);
        out(:,ii) = out(:,ii) + Gauss_noise;
        %---------------------------------------------------------------

             
     end
     

    % reset the global model to update it from the new local models
    out_central=zeros(num_feature,1);  

    % update the global model parameters (capital Theta)
    for ii =1:no_workers
        out_central = out_central +  1/(no_workers)*(out(:,ii)+lambda(:,ii)/rho);
    end
     
         
    % update the dual variables as in Eq. (12)
    for ii=1:no_workers

      lambda(:,ii) = lambda(:,ii) + rho*(out(:,ii)-out_central);

    end
    
    
         
        final_obj = 0.5*norm(XX * out_central - YY)^2;
        

        obj_DPQADMM(i)=final_obj;
        loss_DPQADMM(i)=abs(final_obj-obj0);     % obj0 is the optimal global model obtained by running the fn opt_sol_closedForm       
        
       
        
        % min_loss = loss_DPQADMM(i)

       
        

 end % end for i =1: maxiter   
    

     clc
     msg = ['Simulation finished for "ADMM+DP", moving forward ... '];
     disp(msg)
end % end function
     




