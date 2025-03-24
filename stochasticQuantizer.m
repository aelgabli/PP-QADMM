function [quantized,number_of_bits_toSend]=stochasticQuantizer...
                                            (quantized,current,prev,bitsToSend)


b=bitsToSend;

tau=1/(2^b-1); % this is 1/ number of quantization levels
number_of_bits_toSend =32+length(current)*b; % the number of bits to send the value of R.

diff = current - prev;
% prev
R = max(abs(diff)); % the L-inifinity norm gives the max of the absolute of the vector components
   
% Stochastic Quantization
% Q is eq. (8) in the paper (diff+R)/Delta
% the denominator 2R*tau = 2R/num of quant. levels, which is the quant. step size (Delta)
Delta = 2*R*tau;
Q = (diff+R)/Delta; % R is added to ensure nonnegativity of the quantized values.

p = (Q-floor(Q)); % the quant. prob. as in eq. (12)

% below is the implementation of eq. (9) for stochastic quantization
for i = 1:length(current)
    temp=rand;
    if(temp <=p(i))
        Q(i)=ceil(Q(i));
    else
        Q(i)=floor(Q(i));
    end    
end

quantized = quantized + Delta*Q - R; % implementation  of eq. (15)

end % end function    

