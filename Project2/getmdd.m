function mdd = getmdd(prices)
%%%% Calculate the maximum drawdown from the net value record

n = length(prices);
max_price = prices(1);
mdd = 0;
for i = 1:n
    max_price = max(max_price, prices(i));
    mdd = max(mdd, (max_price - prices(i)) / max_price);
end