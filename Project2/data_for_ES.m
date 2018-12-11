function [mu,P_S_R,w] = data_for_ES(Prices,trade_date, ...
    horizon,sample_frequency,number_of_samples,rate_of_decay)
    
% As is the case when risk is measured using a covariance matrix, 
% for expected shortfall there need to be considerably more samples
% than there are assets in the portfolio, else the solver can conclude
% that some choices of portfolios have little risk but lots of potential
% for return. Consquently, when working with expected shortfall, it is
% useful to have a file similar to adapted_stats.m, where shorter
% frequency samples are modified to act as if they are longer frequency samples.
% This is such a file. It is a slight modfication of adapted_stats.m
% (indeed, only the last line has been added).
% Note that mu is still an output -- indeed, expected shortfall is
% about risk, not about return. Also, the probability measure w for the samples
% is now output, as is a new matrix P_S_R (explained at the end of the file).

% (The many comments appearing in adapted_stats are removed here.)
    
%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = horizon; P = Prices; t_d = trade_date; s_f = sample_frequency;
n_s = number_of_samples; r_d = rate_of_decay;

s_d = t_d - s_f*[0:n_s];    % reverse chronological
s_d = fliplr(s_d);              % chronological
% sample dates, a row vector

S_P = P(s_d,:);
% Sample Prices, a matrix

S_C_R = log(S_P(2:end,:)./S_P(1:(end-1),:));
% Sample Compound Returns, a matrix

% now it ’s time to construct the weight
w=(1-r_d).^(1:n_s);
w=fliplr(w)'/sum(w);
% weights, a non-negative vector that sums to 1

mean_c_r = S_C_R'*w;
% mean vector of compound returns

Cov_C_R = (S_C_R'*diag(w)*S_C_R) - mean_c_r*mean_c_r';
% covariance matrix of compound returns

adapted_mean_c_r = (h/s_f)*mean_c_r;
% adapting mean vector to reflect length of holding period

Adapted_Cov_C_R = (h/s_f)*Cov_C_R;
% adapting covariance matrix to length of holding period

muu = exp(adapted_mean_c_r + .5*diag(Adapted_Cov_C_R)); 
mu = muu - 1;
% resulting mean vector of ARITHMETIC returns

V = (muu*muu').*(exp(Adapted_Cov_C_R) - 1);
% resulting covariance matrix of ARITHMETIC returns

P_S_R = exp((h/s_f)*S_C_R) - 1;
% A synthetic array of arithmetic sample returns obtained by appropriately
% modifying arithmetic sample returns obtained for short time periods
% so as to be useful in a model in which the portfolio will be rebalanced
% only at longer intervals of time.

