close all
clear
clc
%%
% ce,ck1, ck2, phie : parameters, conditioned at age 45
% result = AAA_main(0.01, 0.01, run_day, 'test','true');

run_day = 15;

para_info(1).mean = 88.38;
para_info(1).std = 51.1; % 51.1
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
coll_pts_ = xi_mat(1,:);

%%                    Create set of collocation points
% Condition on: no points to be negative
shift_count = 1;
while (count<1+poly_order*np)
    for i=1:np
        temp = para_value;
        temp(:,i) = circshift(temp(:,i),-shift_count);
        coll_pts = [coll_pts;temp(1,:)];
        
        temp = xi_mat;
        temp(:,i) = circshift(temp(:,i),-shift_count);
        coll_pts_ = [coll_pts_;temp(1,:)];
        
        count = count + 1;
    end
    shift_count = shift_count + 1;
end

for i=1:size(coll_pts,1)
    % --- RHS
    out = AAA_main(0.01, 0.01, run_day, 'test','true',coll_pts(i,:));
    y_fem(i,1) = max(out.max_diameter);
    
    % --- LHS
    H_temp = 1;
    
    for j=1:np
        for k=1:poly_order
            Her_func = Hermite_poly(k);
            H_temp = [H_temp,Her_func(coll_pts_(i,j))];
        end
    end
    K(i,:) = H_temp;
    fprintf('Create CM equations... %d%%\n',round(i/size(coll_pts,1)*100));
end

nancount = 0;
mean_y = nanmean(y_fem);
for i=1:size(y_fem,1)
    if (isnan(y_fem(i,1))==1)
        y_fem(i,1) = mean_y;
        nancount = nancount + 1;
    end
end
fprintf('There are %d NaNs in y_fem\n',nancount);

CM_coeff = K\y_fem;

%% Run two methods through a span of parameters
para_span = parameter_generate(para_info);
for i=1:size(para_span,1)
    out = AAA_main(0.01, 0.01, run_day, 'test','true',para_span(i,:));
    y_FEM(i,1) = max(out.max_diameter);
    
    y_temp = CM_coeff(1);
    for j=1:np
        for k=1:poly_order
            Her_func = Hermite_poly(k);
            y_temp = y_temp + ...
                CM_coeff((j-1)*poly_order+k+1)*Her_func((para_span(i,j)-para_info(j).mean)/para_info(j).std);
        end
    end
end

%% Run MC on FEM model


%% Run MC on CM model