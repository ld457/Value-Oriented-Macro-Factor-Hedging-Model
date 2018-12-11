function [x0,x] = cvx_ES_hedge(mu0,mu,P_S_R,w,alpha,beta,xx,trans_cost, min_x0, B, exposure, crisis_flag)
% The optimization problem is to maximize portfolio return subject to 
% constraints of expected shortfall, macro hedging and other constraints

n = length(mu);
n_s = length(w);
m = length(alpha);

cvx_begin quiet
	variables x0 x(n) y(n) loss(n_s) total_trans_cost ell
	
	maximize( mu0*x0 + mu'*x )
	
	subject to
    
        loss == -mu0*x0 - P_S_R*x
    
        for i=1:m
            ell+(1/(1-alpha(i)))*w'*max([loss-ell,zeros(n_s,1)],[],2) <= beta(i)
        end
	
        x0 + sum(x) + total_trans_cost == 1
        		
		x == xx + y
		
		trans_cost*sum(abs(y)) <= total_trans_cost
		
        % minimum weight of risk-free asset
		x0 >= min_x0
        
        abs(x - 1/(n+1)) <= .05
        
        % Macroeconomic Hedging
        abs(B * x) <= abs(exposure)
        
        % If there is no crisis warning, we incorporate short sales limit
        if crisis_flag == 0
            x >= 0
        end

        
			
cvx_end	





