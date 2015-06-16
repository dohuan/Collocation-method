function out = parameter_generate(para_info)
    iter = 800;
    spanSize = 10;
    np = size(para_info,2);
    for i=1:np
        para_rand = ...
            para_info(i).mean + para_info(i).std.*randn(iter,1);
        out(:,i) = linspace(min(para_rand),max(para_rand),spanSize)';
    end
end