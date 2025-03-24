function [w_optimal, out] = opt_sol_closedForm(XX,YY)
w_optimal = (XX'*XX)\(XX'*YY);

out =0.5*norm(XX*w_optimal - YY)^2;
 




