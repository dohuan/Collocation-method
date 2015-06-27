function out = ortho_poly(order,dinfo)
%%     Return orthogonal polynomial for arbitrary normal distribution
% dist_type.name dist_type.mean dist_type.std
syms x

for i=1:order
    if (i==1)
        H{i} = matlabFunction(x-ortho_int(@(x)x,@(x)x.*0+1,dinfo));
        h{i} = matlabFunction(H{i}(x)/sqrt(ortho_int(@(x)H{i}(x),@(x)H{i}(x),dinfo)));
    elseif(i==2)
        temp = matlabFunction(x.*h{1}(x));
        H{i} = matlabFunction(x.*h{1}(x)-ortho_int(@(x)temp(x),@(x)h{1}(x),dinfo)*h{1}(x)-...
                                sqrt(ortho_int(@(x)H{1}(x),@(x)H{1}(x),dinfo)));
        h{i} = matlabFunction(H{i}(x)/sqrt(ortho_int(@(x)H{i}(x),@(x)H{i}(x),dinfo)));
    else
        temp = matlabFunction(x.*h{i-1}(x));
        H{i} = matlabFunction(x.*h{i-1}(x)-ortho_int(@(x)temp(x),@(x)h{i-1}(x),dinfo)*h{i-1}(x)-...
                       sqrt(ortho_int(@(x)H{i-1}(x),@(x)H{i-1}(x),dinfo))*h{i-2}(x));
        h{i} = matlabFunction(H{i}(x)/sqrt(ortho_int(@(x)H{i}(x),@(x)H{i}(x),dinfo)));
    end
    out.H = H;
    out.h = h;
end

end

function out = ortho_int(p1,p2,dinfo)
    syms x sig mu
    d_func = @(x)1./(sig*sqrt(2*pi))*exp(-((x-mu).^2)/(2*sig^2));
    lim = [-inf inf];
    
    temp = matlabFunction(p1(x)*p2(x)*d_func(x));
    out = integral(@(x)temp(dinfo.mean,dinfo.std,x),lim(1),lim(2));
end