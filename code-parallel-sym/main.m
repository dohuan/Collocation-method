close all
clear
clc
%%
% ce,ck1, ck2, phie : parameters, conditioned at age 45
% result = AAA_main(0.01, 0.01, run_day, 'test','true');

run_day = 100;

para_info(1).mean = 88.38;
para_info(1).std = 51.1; % 51.1
para_info(1).name = 'ce';
para_info(1).dist_type = 'Normal';
para_info(1).trunc = [10 200];

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

% para_info(1).mean = 0.203;
% para_info(1).std = 8.2e-2;
% para_info(1).name = 'phie';
% para_info(1).dist_type = 'Normal';

np = size(para_info,2);
poly_order = 3;

for i=1:np
    %[Ortho{i},para_value(:,i)] = coll_points_generate(poly_order+1,para_info(i));
    [Ortho{i},para_value(:,i)] = coll_points_generate(poly_order,para_info(i));
end
%colPtsNeeded = np*poly_order+1;
%nxi = size(para_value,1);

%%                    Create set of collocation points
% --- Condition on: no points to be negative
% --- Delete all negative value in para_value

% count = 1;
% coll_pts = [];
% shift_count = 1;
% offset = 0;
% break_flag = 0;
% while(break_flag == 0)
%     for i=1:np
%         temp = para_value;
%         for j=1:nxi
%             if (count~=1)
%                 % --- Shift one unit
%                 temp(:,i) = circshift(temp(:,i),-j);
%             end
%             if (isempty(find(temp(1+offset,:)<0, 1))==1)
%                 coll_pts = [coll_pts;temp(1+offset,:)];
%                 count = count + 1;
%             end
%         end
%     end
%     offset = offset + 1;
%     if (offset>nxi-1)
%         fprintf('Generate all possible collocation points:\n');
%         fprintf('Need: %d \n',colPtsNeeded);
%         fprintf('Have: %d \n',size(unique(coll_pts,'rows'),1));
%         break_flag = 1;
%     end
% end
% if (colPtsNeeded>size(unique(coll_pts,'rows'),1))
%     error('Not enough collocation points!');
% end
% coll_pts = unique(coll_pts,'rows','stable');
% coll_pts = coll_pts(1:colPtsNeeded,:); 

coll_pts = [];

for k1=1:poly_order+1
    for k2=1:poly_order+1
        for k3=1:poly_order+1
            for k4=1:poly_order+1
                coll_pts = [coll_pts;[para_value(k1,1) para_value(k2,2)...
                                       para_value(k3,3) para_value(k4,4)]];
            end
        end
    end
end

% --- Create worker pool (comment for local machine)
% N = 5;
% poolobj = gcp('nocreate');
% if isempty(poolobj)
% 	poolsize = 0;
% else
% 	poolsize = poolobj.NumWorkers;
% end
% 
% if poolsize == 0
% 	parpool('local',N);
% else
% 	if poolsize~=N
% 		delete(poolobj);
% 		parpool('local',N);
% 	end
% end

y_fem = zeros(size(coll_pts,1),1);
parfor i=1:size(coll_pts,1)
    % --- RHS
    out = AAA_main(0.01, 0.01, run_day, 'test','false',coll_pts(i,:));
    y_fem(i,1) = max(out.max_diameter);
    
    % --- LHS
    H_temp = [];
    
    for k1=1:poly_order+1
        for k2=1:poly_order+1
            for k3=1:poly_order+1
                for k4=1:poly_order+1
                    H_temp = [H_temp,Ortho{1}.H{k1}(coll_pts(i,1))*Ortho{2}.H{k2}(coll_pts(i,2))*...
                                            Ortho{3}.H{k3}(coll_pts(i,3))*Ortho{4}.H{k4}(coll_pts(i,4))];
                end
            end
        end
    end
    
%     for j=1:np
%         for k=1:poly_order
%             %ortho_func = Hermite_poly(k);
%             
%             ortho_func = Ortho{j}.H{k};
%             H_temp = [H_temp,ortho_func(coll_pts(i,j))];
%             %H_temp = [H_temp,Her_func(coll_pts(i,j))];
%         end
%         % --- add cross-product term
% %         if (j<np)
% %             Her_func = Hermite_poly(1);
% %             H_temp = [H_temp,...
% %                       Her_func(coll_pts_(i,j))*Her_func(coll_pts_(i,j+1))];
% %         else
% %             Her_func = Hermite_poly(1);
% %             H_temp = [H_temp,...
% %                       Her_func(coll_pts_(i,j))*Her_func(coll_pts_(i,1))];
% %         end
%     end
    K(i,:) = H_temp;
    fprintf('Create CM equations... %d%%\n',round(i/size(coll_pts,1)*100));
end
delete(poolobj);

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
% para_span = parameter_generate(para_info);
% for i=1:size(para_span,1)
%     out = AAA_main(0.01, 0.01, run_day, 'test','true',para_span(i,:));
%     y_FEM(i,1) = max(out.max_diameter);
%     
%     y_temp = CM_coeff(1);
%     for j=1:np
%         for k=1:poly_order
%             ortho_func = Hermite_poly(k);
%             y_temp = y_temp + ...
%                 CM_coeff((j-1)*poly_order+k+1)*ortho_func((para_span(i,j)-para_info(j).mean)/para_info(j).std);
%         end
%     end
% end

% --- Generate parameters that have inverse Gaussian distribution
iter = 100; % MCMC runs 
for i=1:np
    mu_ = para_info(i).mean;
    std_ = para_info(i).std;
    pd =  makedist('InverseGaussian','mu',mu_,'lambda',mu_^3/std_2);
    %para_test(:,i) = gen_ig(iter,mu_,mu_^3/std_^2); % inverse Gaussian
    para_test(:,i) = random(pd, iter, 1); % inverse Gaussian
end
%% Run MC on FEM model
tic
y_fem_test = zeros(iter,1);
% --- Create worker pool
poolobj = gcp('nocreate');
if isempty(poolobj)
	poolsize = 0;
else
	poolsize = poolobj.NumWorkers;
end

if poolsize == 0
	parpool('local',N);
else
	if poolsize~=N
		delete(poolobj);
		parpool('local',N);
	end
end

for i=1:iter
    out = AAA_main(0.01, 0.01, run_day, 'test','false',para_test(i,:));
    y_fem_test(i,1) = max(out.max_diameter);
    fprintf('Run FEM test... %d%%\n',round(i/iter*100));
end
delete(poolobj);

time_FEM = toc/60;
fprintf('Collapsed time for FEM: %.2f mins \n',time_FEM);

%% Run PCM on CM model
tic
for i=1:iter
    H_temp = [];
    for k1=1:poly_order+1
        for k2=1:poly_order+1
            for k3=1:poly_order+1
                for k4=1:poly_order+1
                    H_temp = [H_temp,Ortho{1}.H{k1}(para_test(i,1))*Ortho{2}.H{k2}(para_test(i,2))*...
                        Ortho{3}.H{k3}(para_test(i,3))*Ortho{4}.H{k4}(para_test(i,4))];
                end
            end
        end
    end
    y_cm_test(i,1) = H_temp*CM_coeff;
end
time_CM = toc/60;
fprintf('Collapsed time for CM: %.2f mins \n',time_CM);

save('./results/101615');

%plot(y_fem_test,'b','LineWIdth',2);
%hold on
%plot(y_cm_test,'r','LineWIdth',2);
%hold off
%legend('FEM','CM');
%set(gca,'FontSize',16);

%figure(2);
%err_test = sqrt((y_cm_test-y_fem_test).^2);
%plot(err_test,'b','LineWidth',2);