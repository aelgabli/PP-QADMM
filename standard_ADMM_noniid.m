function [obj_ADMM, loss_ADMM] = standard_ADMM_noniid(XX,YY, rho, no_workers, ...
        num_feature, noSamples, num_iter, obj0)

Iter= num_iter;              
lambda = zeros(num_feature,no_workers);
out=zeros(num_feature,no_workers);
out_central=zeros(num_feature,1);
%prev_out=zeros(s1,no_workers);
%q_out=zeros(s1,no_workers);
% snr=10^(20.5/10);

max_iter = num_iter;
% energy(1)= 0;

     clc
     msg = ['Running Simulation for "ADMM", hold on ... (1 out of 5) '];
     disp(msg)


 for i = 1:max_iter
     % if(i > 1)
     %     energy(i)=energy(i-1);
     % end

     clc
     msg = ['Itaeration # ', num2str(i), ' for "ADMM", hold on ... (1 out of 5) '];
     disp(msg)
         
     
     for ii =1:no_workers

         term_1=rho*out_central;

         B1 = lambda(:,ii);
         
         first = sum(noSamples(1:ii-1))+1;
         last = first+noSamples(ii)-1;
        
        X=XX(first:last,1:num_feature);
        Y=YY(first:last);

        % The following equation calculates the local model based on the
        % average loss functions, not their sum.
        x=((1/noSamples(ii))*(X'*X)+rho*eye(num_feature))\((1/noSamples(ii))*(X'*Y)-B1+term_1);
        %x=(X'*X+rho*eye(num_feature))\(X'*Y-B1+term_1); % this is based on the sum of loss functions
        
        out(:,ii) =x;
     
                
        
     end
     

    
    out_central=zeros(num_feature,1);
    for ii =1:no_workers
       out_central = out_central + 1/(no_workers)*out(:,ii); %(out(:,ii)+lambda(k,ii)/rho);
    end
    
         
    
    for ii=1:no_workers

      lambda(:,ii) = lambda(:,ii) + rho*(out(:,ii)-out_central);

    end
    
        final_obj=0;
        for j=1:no_workers
            first = sum(noSamples(1:j-1))+1;
            last = first+noSamples(j)-1;
                
            X=XX(first:last,1:num_feature);
            Y=YY(first:last);
            final_obj =final_obj+1/noSamples(j)*0.5*norm(X * out_central - Y)^2;% we already compute the average loss
        end
    
         
       
        %final_obj = 0.5*norm(XX * out_central - YY)^2;
        
        obj_ADMM(i)=final_obj;
        loss_ADMM(i)=abs(final_obj-obj0);            
        
        
       
    end   
    
     clc
     msg = ['Simulation finished for "ADMM", moving forward ... '];
     disp(msg)

end
     




