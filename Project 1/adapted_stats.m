function [mu,V] = adapted_stats(Prices,trade_date, ...
    horizon,sample_frequency,number_of_samples,rate_of_decay)

%mu: mean vector of historical returns

%V: covariance matrix of historical returns

%Prices: Matrix of asset prices, with rows corresponding to time periods
% (ordered chronologically), and with columns corresponding to assets.

%trade_date: The number of the period prior to which data can be
% accessed for calculating mu and V.  Think of a new portfolio as being
% constructed at the beginning of the period whose number is trade_date.
% The only prices then known would be ones from earlier periods.)

%horizon: The number of periods which the new portfolio will be held.
% (For example, if the portfolio is to be held one month and the data
% is weekly, then horizon=4, whereas if data is daily, horizon=20.)

%sample_frequency: The number of periods between sampled returns in 
% calculating mu and V.  For example, if the historical returns are
% to be sampled weekly, then sample_frequency=1 if data is weekly,
% whereas sample_frequency=5 if data is daily.

%number_of_samples: The number of historical returns to compute. Must
% satisfy the condition 
%          sample_frequency*number_of_samples < trade_date  
% since otherwise there will insufficient historical data.

%rate_of_decay: Each of the historical returns will be assigned a weight.
% The most recent return will be given the greatest weight.  The
% weight of the next most recent return will be (1-rate_of_decay) times
% the weight of the most recent return.  Likewise, the weight of the
% second most recent return will be (1-rate_of_decay)^2 times the
% weight of the most recent return.  And so on.
% The weights will sum to 1, thus giving a probability measure on the 
% samples.  This is the probability measure used in computing mu and V.
% Be careful to satisfy the condition 0 <= rate_of_decay < 1.
% (Clearly, if all samples are to be given the same weight, choose
% rate_of_decay=0.
% BE CAREFUL to change rate_of_decay when changing sample_frequency.
% In particular, if r1 is the rate chosen for frequency f1, then for
% frequency f2 the equivalent rate is 
%                  r2 = 1 - (1-r1)^(f2/f1)
% These two rates would result in giving the same relative weight 
% to a sample return from one month in the past, one year in the past, 
% or however long ago. In other words, if r1 causes samples from
% one year ago to be given half as much weight as the most recent sample,
% so will r2.

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

