function out = ortho_poly(order,dinfo)
%%     Return orthogonal polynomial for arbitrary normal distribution
% dist_type.name dist_type.mean dist_type.std
syms x

if (order==0)
    out.H = 1;
    out.h = 1;
else
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

end

function out = ortho_int(p1,p2,dinfo)
syms x   % sig mu
mean = dinfo.mean;
sig = dinfo.std;
dist_type = 1; % 0: normal dist 1: inverse Gamma
if(dist_type==0)
    d_func = @(x)1./(sig*sqrt(2*pi))*exp(-((x-mean).^2)/(2*sig^2));
    lim = [-inf inf];
elseif(dist_type==1)
    %alpha = mean^2/sig+2;
    %beta = mean*(mean^2/sig+1);
    %d_func = @(x)(mean.*(mean.^2./sig+1)).^(mean.^2./sig+2)./...
    % gamma(mean.^2./sig+2).*x.^(-(mean.^2./sig+2)-1).*exp(-(mean.*(mean.^2./sig+1)).*x.^(-1));
    d_func = @(x)pdf('InverseGaussian',x,mean,mean^3/sig^2);
    lim = [0 inf];
end
temp = @(x)p1(x).*p2(x).*d_func(x);
%temp = matlabFunction(p1(x)*p2(x)*d_func(x));
%out = integral(@(x)temp(dinfo.mean,dinfo.std,x),lim(1),lim(2));
out = integral(@(x)temp(x),lim(1),lim(2));
end