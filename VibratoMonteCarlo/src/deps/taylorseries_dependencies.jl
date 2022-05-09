using .TaylorSeries

v_value(x::AbstractSeries) = x[0]

function bell_n(X)
	if isempty(X)
		return 1.0;
	end
	result=0.0
	n=length(X)-1
	for i=0:n
		coeff_bin_i=binomial(n,i);
		@views X_i=X[1:(n-i)];
		@views X_i_1=X[i+1]
		prev_bell=bell_n(X_i)
		result+=coeff_bin_i*prev_bell*X_i_1
	end
	return result
end

function faa_di_bruno_bell(x::Taylor1);
	max_ord=get_order(x)
	orders=1:max_ord;
	@views derivatives=[differentiate(x,i)[0] for i in orders]
	@views vec=[bell_n(derivatives[1:i])/factorial(i) for i in orders]
	# @views vec=[bell_n(derivatives[1:i]) for i in orders]
	# check_2=derivatives[1]^2+derivatives[2]
	# @show check_2-vec[2]
	# check3_=derivatives[3]+3*derivatives[1]*derivatives[2]+derivatives[1]^3
	# @show check3_-vec[3]*factorial(2)
	return Taylor1([1.0,vec...]);
	# return [derivatives[1],check_2,check3_];
	# return vec
end

function v_mod(x::AbstractSeries)
	return faa_di_bruno_bell(x);
end