% function  multi_period_with_ES

load data.mat; 

% 19 stocks that can be explained by at least one macro factor
% Price = Price(:, [1,3,4,7,8,10,11,12,13,14,15,16,17,18,19,20,22,24,26]);

[m,n] = size(Price);
	% n = number of risky assets
e = ones(n,1);	

% Consumer sentiment, warning system indicator
sentiment = factor((start-2):end,8);

% 5 factors that have high explanatory power
% factor = factor(:, [1, 3, 4, 8, 9]); 
[~, factor_num] = size(factor);

%%%% factor
%%% 1 - Average weekly hours
%%% 2 - Average weekly jobless claims for unemployment insurance
%%% 3 - Manufacturers' new orders for consumer goods/materials
%%% 4 - Manufacturers' new orders for non-defense capital goods
%%% 5 - Building permits
%%% 6 - Money Supply (M2)
%%% 7 - Interest rate spread (10-year Treasury vs. Federal Funds target)
%%% 8 - Index of consumer expectations
%%% 9 - S&P 500

%%%%%%%%%%%%%%%
%%%% PARAMETERS

horizon = 1; % rebalance monthly (every months)
start = 4*12; % the month in which you are first given a portfolio to rebalance 

number_of_samples = start - 1; % how many samples are to be used 
						% in computing return avereages and covariances
                  
number_rebalances = m - number_of_samples - 1; % the number of times the portfolio will be rebalanced
sample_frequency = 1; % 1 = monthly.
						
r_w_f_o_y_e = .6; % "relative weight for one year earlier" 
				   % -- a value .4 means that for the (exponential) weights 
				    % used in computing return averages and covariances, 
				     % the weight assigned to the time period one year ago
				      % should be .4 times the weight assigned 
				       % to the most recent period.   

alpha = [.99,.95,.85];
% Each entry is a confidence level at which expected shortfall is to be
% controlled


allowable_risk = [1, 1.25, 1.5];
% allowable_risk = [2, 2.5, 3];
% allowable_risk = [0.5,0.625,0.75];
% allowable_risk = [5, 6.25, 7.5];
% allowable_risk = [0.2,0.25,0.3];

% allowable_risk = [1.5, 1, 1.25];
% allowable_risk = [1.25, 1.5, 1];
    % For the various confidence levels, these are the amounts of expected
    % shortfall you want not to exceed, where the measures are relative to
    % the expected shortfall of the benchmark portfolio. For example, the
    % first entry being 1 means that for the most extreme .01% of tail events, risk
    % no greater than the benchmark risk is to be assumed, whereas the 1.5 in
    % the last entry means that for the tail consisting of the 15% worst
    % scenarios (most of these being much milder than the extreme .01% of
    % events above), we are willing to take on 1.5 times the benchmark risk.

allowable_exposure = 1; % This is the level of exposure to macro factors relative to the benchmark portfolio

%%%%%% Crisis risk management parameters

% The allocation weight of risk-free asset when consumer sentiment touches red line
crisis_x0 = 0.8;
% The allocation weight of risk-free asset in normal economy
normal_x0 = 0.0;
% crisis threshold
crisis_threshold = -0.04;
% recovery threhold
recov_threshold = 0.03;
% crisis flag
crisis_flag = 0;
 
 %%%%%% END OF crisis risk management parameters

trans_cost = .005;  % transaction cost
	
wealth = 10000; % initial wealth measured in dollars, including money invested in assets
			   % (one dollar invested in an asset is considered as one dollar of wealth,
			   %  even though in liquidating the asset, transaction costs would be paid)	
			   
x0 = 1/(n+1); % proportion of wealth in bank initially
x = ones(n,1)/(n+1); % proportions in risky assets initially

% Assume the benchmark portfolio is initally equal-weighted, with 1/(n+1) being the 
% proportion of wealth invested in each asset and in the bank.
									
%%%% END OF PARAMETERS	
%%%%%%%%%%%%%%%%%%%%%%

rate_of_decay = 1 - r_w_f_o_y_e^(sample_frequency/52);

initial_wealth = wealth;
benchmark_wealth = wealth;	

rebalance_dates = start + horizon*(0:number_rebalances-1);

my_wealth_record = zeros(length(rebalance_dates)+1, 1);
my_wealth_record(1) = initial_wealth;
benchmark_wealth_record = zeros(length(rebalance_dates)+1, 1);
benchmark_wealth_record(1) = benchmark_wealth;

rf_record = zeros(length(rebalance_dates), 1);
rf_position = zeros(length(rebalance_dates), 1);

flags = zeros(length(rebalance_dates)+1, 1);
flags(1) = crisis_flag;

for i = 1:length(rebalance_dates)

	trade_date = rebalance_dates(i);

    %%%%% REGRESSION ON MACRO FACTOR TO GET THE EXPOSURES
    
    s_d = trade_date - sample_frequency*[0:number_of_samples];    % reverse chronological
    s_d = fliplr(s_d);              % chronological
    % sample dates, a row vector

    S_P = Price(s_d,:);
    % Sample Prices, a matrix

    S_R = S_P(2:end,:)./S_P(1:(end-1),:) - 1; % Sample return, a matrix
    
    % X = [ones(number_of_samples, 1), factor(trade_date - number_of_samples+1:trade_date, :)];
    X = factor(trade_date - number_of_samples+1:trade_date, :);
    B = zeros(factor_num, n);
    
    for stock = 1:n
        % Regression
        lm = fitlm(X, S_R(:, stock));
        pvalues = lm.Coefficients.pValue(2:end);
        coef = lm.Coefficients.Estimate(2:end);
        
        % Limit factors to be the ones with relatively high explanatory
        % power
        coef(pvalues > 0.1) = 0; 
%         [~, ind] = sort(pvalues, 'descend');
%         coef(ind(1:4)) = 0;

        B(:, stock) = coef;
    end
    
    benchmark_beta = sum(B, 2) / (n+1);
    exposure = benchmark_beta * allowable_exposure;  
    
    %%%%% END OF REGRESSION ON MACRO FACTOR TO GET THE EXPOSURES
    
    %%%%% REBALANCE YOUR PORTFOLIO AND PAY TRANSACTION COSTS %%%%%%
    % It is more natural to rebalance the benchmark portfolio later %
		
	[mu,P_S_R,w] = data_for_ES(Price,trade_date, ...
    	horizon,sample_frequency,number_of_samples,rate_of_decay);
    
    mu0 = (1+.01*risk_free_rate(trade_date-1))^(horizon/52) - 1;
    rf_record(i) = mu0; % record the risk free rate in every period
    
    benchmark_risk = calculate_expected_shortfall(mu0,P_S_R,w,1/(n+1),ones(n,1)/(n+1),alpha);
    								     % including the bank		
    beta = allowable_risk.*benchmark_risk; 
    
    xx = x;    
    
    %%%%%%% Sentiment warning risk management
    if sentiment(i+1) < crisis_threshold && sentiment(i) < crisis_threshold && sentiment(i+2) < crisis_threshold
        crisis_flag = 1;
        [x0,x] = cvx_ES_hedge(mu0,mu,P_S_R,w,alpha,beta,xx,trans_cost,crisis_x0, B,exposure, crisis_flag);
    elseif sentiment(i+1) > recov_threshold && sentiment(i+2) > recov_threshold
        crisis_flag = 0;
        [x0,x] = cvx_ES_hedge(mu0,mu,P_S_R,w,alpha,beta,xx,trans_cost,normal_x0, B,exposure, crisis_flag);
    elseif crisis_flag == 1
        [x0,x] = cvx_ES_hedge(mu0,mu,P_S_R,w,alpha,beta,xx,trans_cost,crisis_x0, B,exposure, crisis_flag);
    else
        [x0,x] = cvx_ES_hedge(mu0,mu,P_S_R,w,alpha,beta,xx,trans_cost,normal_x0, B,exposure, crisis_flag);
    end
    flags(i+1) = crisis_flag;
    
    %%%%%%% End of sentiment warning risk management
    
%     [x0,x] = cvx_ES_hedge(mu0,mu,P_S_R,w,alpha,beta,xx,trans_cost,normal_x0, B,exposure, 1);
    rf_position(i) = x0;
    
    wealth = wealth*(x0 + sum(x));
    	% This is the same thing as updating your wealth by subtracting
    	% all transaction costs from the rebalancing.  Indeed, in rebalancing,
    	% the proportion of your wealth going to trans costs is 1 - x0 - sum(x).
        
    total = x0 + sum(x);
    x0 = x0/total;
    x = x/total;
    	% Rescaling x0 and x so that the sum is 1 (i.e., proportions of current wealth)
    	
    %%%%%% PROCEED TO END OF TIME PERIOD AND ACCOUNT FOR GAINS, LOSSES %%%%%	
    	  
    returns = (Price(trade_date+horizon,:)-Price(trade_date,:))./Price(trade_date,:);
    	% vector of actual returns for risky assets (this is a row vector)	
    
    multiplier = 1 + mu0*x0 + returns*x;	
    wealth = multiplier*wealth;
    my_wealth_record(i+1) = wealth;
        	% by leaving off the semicolon, you can watch how wealth changes as the program runs

    if wealth<=0
     break; % stops the program if bankruptcy occurs
     		% Not needed for benchmark portfolio (because it is long only)
    end
    
    x0 = (1+mu0)*x0/multiplier;
    x = x.*(1+returns)'/multiplier;
    % these are the proportions of current wealth invested in assets
    
    % Now its time to rebalance the benchmark portfolio and pay transaction costs
    
    benchmark_x0 = (1+mu0)/(n+1);
    benchmark_x = (1+returns)/(n+1);
    % This gives how the equal-weighted portfolio has changed during the time period.
    % The initial unit of wealth has become  benchmark_x0 + sum(benchmark_x).
    % This new level of wealth needs to be distributed equally among the assets and bank.
    % The optimal amount z to put into each one is determined by the following function,
    % which finds the value z so as to minimize transaction costs
    
    z = rebalance_benchmark(benchmark_x0,benchmark_x,trans_cost);  	
	benchmark_wealth = benchmark_wealth*(n+1)*z;
    benchmark_wealth_record(i+1) = benchmark_wealth;
		
	% Until the end of the next time period, 
	% think of the benchmark portfolio as having been rebalanced
	% with wealth divided equally, that is, the portion of wealth invested in
	% each asset and the bank is 1/(n+1).	
    
    % your_wealth_vs_benchmark = [wealth, benchmark_wealth]
    
   
    
end

%%%%% Performance of Strategy
my_return = my_wealth_record(2:end) ./ my_wealth_record(1:end-1) - 1;
bchmk_return = benchmark_wealth_record(2:end) ./ benchmark_wealth_record(1:end-1) - 1;

% annual return
my_annual_return = (wealth/initial_wealth)^(12/(horizon*number_rebalances))-1;
bchmk_annual_return = (benchmark_wealth/initial_wealth)^(12/(horizon*number_rebalances))-1;

% volatility
my_std = std(my_return) * sqrt(12);
bchmk_std = std(bchmk_return) * sqrt(12);

% sharpe ratio
my_sharpe = mean(my_return - rf_record)*12 / my_std;
bchmk_sharpe = mean(bchmk_return - rf_record)*12 / bchmk_std;

% maximum drawdown
my_mdd = getmdd(my_wealth_record);
bchmk_mdd = getmdd(benchmark_wealth_record);

fprintf('your final bank account %f\n');
x0
fprintf('your final risky portfolio %f\n');
x

fprintf('your final wealth %f\n',wealth);
fprintf('benchmark final wealth %f\n',benchmark_wealth);

fprintf('your annualized rate of return %f\n', my_annual_return);
fprintf('your annualized volatility %f\n', my_std);
fprintf('your sharp ratio %f\n', my_sharpe);
fprintf('your maximum drawdown %f\n',my_mdd);

fprintf('benchmark annualized rate of return %f\n', bchmk_annual_return);
fprintf('benchmark annualized volatility %f\n', bchmk_std);
fprintf('benchmark sharp ratio %f\n', bchmk_sharpe);
fprintf('benchmark maximum drawdown %f\n',bchmk_mdd);

date = datetime(2003,12,1);
%date = datetime(2008,12,1);
dates = date+ calmonths(0:number_rebalances);

sp = sp500(start:start+number_rebalances)./sp500(start) * initial_wealth;
plot(dates, my_wealth_record, dates, benchmark_wealth_record,dates, sp)
title(['Strategy Performance (Smaller Stock Pool with 5 Factors)'])
xlabel('Time')
ylabel('Wealth')
legend('Hedging Model Strategy','Benchmark Strategy', 'SP500','Location','northwest')
