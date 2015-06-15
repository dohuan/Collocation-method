tic;

%% Setup
mu = [5:5:25];
sigma = [1:2:10];
k = [0.5:0.05:0.8];

% %% Constant mu
% mu_const = 15;
% count = 0;
% for sigma_index = 1:size(sigma,2)
%     for k_index = 1:size(k,2)
%         count = count + 1;
%         d.mu = mu_const;
%         d.sigma = sigma(sigma_index);
%         d.k = k(k_index);
%         params(count) = d;
%     end
% end
% parfor index = 1:count
%     damage_params = params(index);
%     namestr = ['muconst_sigma' num2str(damage_params.sigma) '_k' num2str(damage_params.k)];
% %     disp(namestr);
%     AAA_main(0.01, 0.01, 30, namestr, true, damage_params);
% end

%% Full Grid
count = 0;
for mu_index = 1:size(mu,2)
    for sigma_index = 1:size(sigma,2)
        for k_index = 1:size(k,2)
            count = count + 1;
            d.mu = mu(mu_index);
            d.sigma = sigma(sigma_index);
            d.k = k(k_index);
            params(count) = d;
        end
    end
end
parfor index = 1:count
    damage_params = params(index);
    namestr = ['mu' num2str(damage_params.mu) '_sigma' num2str(damage_params.sigma) '_k' num2str(damage_params.k)];
%     disp(namestr);
    AAA_main(0.01, 0.01, 10, namestr, true, damage_params);
    disp(['Finished run ' namestr]);
end


toc;
