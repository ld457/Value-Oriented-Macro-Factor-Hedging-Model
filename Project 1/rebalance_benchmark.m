function z = rebalance_benchmark(benchmark_x0,benchmark_x,trans_cost)

n = length(benchmark_x);

    cvx_begin quiet
        variable z
	
        maximize(z)
        subject to
            (n+1)*z <= benchmark_x0 + sum(benchmark_x) - trans_cost*sum(abs(z - benchmark_x))
cvx_end


