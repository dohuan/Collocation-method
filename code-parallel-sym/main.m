close all
clear
clc
%%
% ce,ck1, ck2, phie : parameters, mean takes at age 45
% run_day = 30;
% result = AAA_main(0.01, 0.01, run_day, 'test','true');
para_info(1).mean = 88.38;
para_info(1).std = 51.1;
para_info(1).name = 'ce';
para_info(1).dist_type = 'Normal';

para_info(2).mean = 1361.12;
para_info(2).std = 699.9;
para_info(2).name = 'ck1';
para_info(2).dist_type = 'Normal';

para_info(3).mean = 12.65;
para_info(3).std = 23.3;
para_info(3).name = 'ck2';
para_info(3).dist_type = 'Normal';

para_info(4).mean = 0.203;
para_info(4).std = 8.2e-2;
para_info(4).name = 'phie';
para_info(4).dist_type = 'Normal';

np = size(para_info,2);
poly_order = 5;
xi_vec = coll_points_generate(poly_order+1);
nxi = size(xi_vec,1);
xi_mat = repmat(xi_vec,1,np);
% --- para_value: [ ce,ck1, ck2, phie ]
for i=1:np
    para_value(:,i) = xi_vec*para_info(i).std + para_info(i).mean;
end
count = 2;
coll_pts = para_value(1,:);
shift_count = 1;
while (count<numel(para_value)+1)
    for i=1:np
        temp = para_value;
        temp(:,i) = circshift(temp(:,i),-shift_count);
        coll_pts = [coll_pts;temp(1,:)];
        count = count + 1;
    end
end




%% Run MC on FEM model


%% Run MC on CM model