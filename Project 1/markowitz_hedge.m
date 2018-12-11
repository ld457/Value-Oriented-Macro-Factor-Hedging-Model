function [x0,x] =  markowitz_hedge(mu0,mu,V,sigma,xx0,xx,trans_cost,B, beta);

n = length(mu);
U = chol(V);


    cvx_begin quiet
    
        variables x0 x(n) y(n) total_trans_cost
              
        maximize(mu0*x0 + mu'*x)
        
        subject to
                    
                norm(U*x) <= sigma
                
                x0 + sum(x) + total_trans_cost == 1
                
                x == xx + y
                
                trans_cost*sum(abs(y)) <= total_trans_cost
                
                x0 >= 0
                
                abs(B * x) <= abs(beta)
                
                abs(x - 1/(n+1)) <= .05
                
                
                
                %.05 <= x <= .2
               
                
               % x(2) + x(5) <= x(6) + x(7) + x(9)
                
                % sum_smallest(x,4) >= -1.25
                    
                    
    cvx_end