%function backtest

load data.mat; 

factor = factor(:, [1,2,3,4,5,7,9,12]);
%factor = factor(:, [1,2,3,4,5,9,12]);
%factor = factor(:, [1,2,3,4,5,7,8,9,10,12]);
%factor = factor(:, [1,2,3,4,5,6,7,9,11,12]);
%factor = factor(:, [1,2,4,5,7]);

%%%% factor
%%% 1 - cpi
%%% 2 - gdp
%%% 3 - federal fund rates
%%% 4 - 10-year Treasury rates
%%% 5 - SP500
%%% 6 - crude oil
%%% 7 - unemployment rate
%%% 8 - industrial production
%%% 9 - dollar index
%%% 10- home price
%%% 11- gold
%%% 12- vix

[m,n] = size(Price);
	% n = number of risky assets
    % m = number of months
e = ones(n,1);	

%%%%%%%%%%%%%%%
%%%% PARAMETERS
%%%% (the following choices of parameters can easily be changed)

horizon = 1; % rebalance monthly (every months)
start = 4*12; % the month in which you are first given a portfolio to rebalance
 
number_of_samples = start - 1; % how many samples are to be used 
						% in computing return avereages and covariances    
                        
number_rebalances = m - number_of_samples - 1; % the number of times the portfolio will be rebalanced 
%number_rebalances = 42;

sample_frequency = 1; % 1 = monthly 

r_w_f_o_y_e = .6; % "relative weight for one year earlier" 
				   % -- a value .4 means that for the (exponential) weights 
				    % used in computing return averages and covariances, 
				     % the weight assigned to the time period one year ago
				      % should be .4 times the weight assigned 
				       % to the most recent period.    	 

allowable_risk = 1;
	% This is the level of risk relative to the benchmark portfolio,
	%   where risk is measured as standard deviation of portfolio returns.
	% Choosing this value to equal 1 means exactly the same amount of risk is allowed,
	% whereas choosing 2 means twice as much risk is allowed as the benchmark, and so on.

allowable_exposure = 1; % This is the level of exposure to macro factors relative to the benchmark portfolio
							
trans_cost = .005;  % transaction cost
	
wealth = 10000; % initial wealth measured in dollars, including money invested in assets
			   % (one dollar invested in an asset is considered as one dollar of wealth,
			   %  even though in liquidating the asset, transaction costs would be paid)	
			   
x0 = .3; % proportion of wealth in bank initially
x = (.7/n)*e; % proportions in risky assets initially

% Assume the benchmark portfolio is initally equal-weighted, with 1/(n+1) being the 
% proportion of wealth invested in each asset and in the bank.
									
%%%% END OF PARAMETERS	
%%%%%%%%%%%%%%%%%%%%%%

rate_of_decay = 1 - r_w_f_o_y_e^(sample_frequency/12);

initial_wealth = wealth;
benchmark_wealth = wealth;	

rebalance_dates = start + horizon*(0:number_rebalances-1);

my_wealth_record = zeros(length(rebalance_dates)+1, 1);
my_wealth_record(1) = initial_wealth;
benchmark_wealth_record = zeros(length(rebalance_dates)+1, 1);
benchmark_wealth_record(1) = benchmark_wealth;

rf_record = zeros(length(rebalance_dates), 1);

for i = 1:length(rebalance_dates)

	trade_date = rebalance_dates(i);
    
    %%%%% REGRESSION ON MACRO FACTOR TO GET THE EXPOSURES
    
    s_d = trade_date - sample_frequency*[0:number_of_samples];    % reverse chronological
    s_d = fliplr(s_d);              % chronological
    % sample dates, a row vector

    S_P = Price(s_d,:);
    % Sample Prices, a matrix

    S_R = S_P(2:end,:)./S_P(1:(end-1),:) - 1; % Sample return, a matrix
    
    X = [ones(number_of_samples, 1), factor(trade_date - number_of_samples+1:trade_date, :)];
    B = X\S_R;
    B = B(2:end, :);
    
    benchmark_beta = sum(B, 2) / (n+1);
    beta = benchmark_beta * allowable_exposure;

    %%%%% REBALANCE YOUR PORTFOLIO AND PAY TRANSACTION COSTS %%%%%%
    % It is more natural to rebalance the benchmark portfolio later %
		
    [mu,V] = adapted_stats(Price,trade_date, ...
    		horizon,sample_frequency,number_of_samples,rate_of_decay);
    %[mu,V] = stats(Price,trade_date,sample_frequency,number_of_samples,rate_of_decay);
    mu0 = (1+.01*risk_free_rate(trade_date))^(horizon/12) - 1;
    rf_record(i) = mu0; % record the risk free rate in every period
    
    benchmark_risk = sqrt(e'*V*e)/(n+1); % there are n+1 financial instruments
    								     % including the bank		
    sigma = allowable_risk*benchmark_risk; 
    
    xx0 = x0;
    xx = x;
    [x0,x] =  markowitz_hedge(mu0,mu,V,sigma,xx0,xx,trans_cost, B, beta);
    
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
title(['Strategy Performance - Full Sample (c = ', num2str(allowable_exposure), ')'])
xlabel('Time')
ylabel('Wealth')
legend('Hedging Model Strategy','Benchmark Strategy', 'SP500','Location','northwest')


