function es = calculate_expected_shortfall(mu0,P_S_R,w,y0,y,alpha)
% Given a FIXED portfolio (y0,y), this function computes the expected
% shortfall at confidence level alpha.

n_s = length(w);
% number of samples

m = length(alpha);
% number of confidence levels for which ES is to be computed

es = zeros(1,m);
for i=1:m
    
    cvx_begin quiet
    
    variables ell loss(n_s)
    
    minimize( ell+(1/(1-alpha(i)))*w'*max([loss - ell, zeros(n_s,1)],[],2) )
    subject to
    loss == -mu0*y0 - P_S_R*y
    cvx_end
    
    es(i) = ell+(1/(1-alpha(i)))*w'*max([loss - ell, zeros(n_s,1)],[],2);
    
end

