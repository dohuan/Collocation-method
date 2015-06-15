function H = Hermite_poly(order)
%% Return the Hermit polynomial
switch (order)
    case 1
        H = @(x)x;
    case 2
        H = @(x)x.^2-1;
    case 3
        H = @(x)x.^3-3*x;
    case 4
        H = @(x)x.^4-6*x.^2+3;
    case 5
        H = @(x)x.^5-10*x.^3+15;
    case 6
        H = @(x)x.^6-15*x.^4+45*x.^2-15;
    case 7
        H = @(x)x.^7-21*x.^5+105*x.^3-105*x;
    case 8
        H = @(x)x.^8-28*x.^6+210*x.^4-420*x.^2+105;
    case 9
        H = @(x)x.^9-36*x.^7+378*x.^5-1260*x.^3+945*x;
end
end