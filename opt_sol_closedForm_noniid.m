function [w_optimal, out] = opt_sol_closedForm_noniid(XX,YY,noSamples,no_workers,num_feature)
temp1=0;
temp2=0;
for i=1:no_workers
    first = sum(noSamples(1:i-1))+1;
    last = first+noSamples(i)-1;
        
    X=XX(first:last,1:num_feature);
    Y=YY(first:last);
    temp1=temp1+(1/noSamples(i))*(X'*X);
    temp2=temp2+(1/noSamples(i))*(X'*Y);
end

w_optimal = (temp1)\(temp2);

out=0;
for i=1:no_workers
    first = sum(noSamples(1:i-1))+1;
    last = first+noSamples(i)-1;
        
    X=XX(first:last,1:num_feature);
    Y=YY(first:last);
    out =out+1/noSamples(i)*0.5*norm(X*w_optimal - Y)^2;
end

     




