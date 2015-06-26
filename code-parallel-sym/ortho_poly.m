function out = ortho_poly(order,dinfo)
%%
% dist_type.name dist_type.mean dist_type.std
syms x
pd = makedist(dinfo.name,dinfo.mean,dinfo.std);
h = @(x)pdf(pd,x);


end

function out = ortho_int(p1,p2,p,lim)
    syms x
    temp = @(x)p1*p2*p;
    temp = matlabFunction(temp);
    out = integral(@(x)temp,[lim(1) lim(2)]);
end