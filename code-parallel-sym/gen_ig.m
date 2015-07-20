function y = gen_ig(n,lambda,mu)

tolerance = 1e-4;

% uniformly distributed random variable
z = rand(1,n);

y = zeros(1,n);

for c = 1:n
    
    lohi = [0 1];
    
    % invert the IG CDF
    % exponential search until something greater than z is found
    
    while (cdf_ig(lohi(2),lambda,mu) < z(c))
        lohi(2) = lohi(2) * 2;
    end
    
    % search the range until found
    while (abs(lohi(2) - lohi(1)) > tolerance)
        mid = (lohi(1) + lohi(2))/2;
        if (cdf_ig(mid,lambda,mu) < z(c))
            lohi(1) = mid;
        else
            lohi(2) = mid;
        end
    end
    
    y(c) = (lohi(1) + lohi(2))/2;
    
end

end

function y = cdf_ig(x,lambda,mu)

y = (1 + erf(sqrt(lambda./(2*x)).*(x/mu - 1)))./2;
y = y + exp(2*lambda/mu)*(1 + erf(-1*sqrt(lambda./(2*x)).*(x/mu + 1)))./2;
end