order = 3;
dinfo.name = 'Normal';
dinfo.mean = 0;
dinfo.std = 1;

out = ortho_poly(order,dinfo)

%--- test orthogonal poly function (compare to Hermite_poly.m
dinfo.name = 'Normal';
dinfo.mean = 0;
dinfo.std = 1;
t = -5:0.1:5;
orth = ortho_poly(9,dinfo);
syms x
for i=1:9
	temp = Hermite_poly(i);
	root_Her = solve(temp(x)==0);
	root_Her = real(double(root_Her));
	
	root_orth = solve(orth.H{i}(x)==0);
	root_orth = real(double(root_orth));
	fprintf('Roots from Hermite: %.2f\n',root_Her);
	fprintf('Roots from Orth: %.2f\n',root_orth);
	%t_Her = temp(t);
	%t_orth = orth.H{i}(t);
	%t_orth = orth.h{i}(t);
	%figure(i)
	%plot(t_Her,'r');
	%hold on
	%plot(t_orth,'b');
	%hold off
	
end

for i=1:size(coll_pts,1)
    % --- LHS
    H_temp = 1;
    for j=1:np
        for k=1:poly_order
            ortho_func = Ortho{j}.H{k};
            H_temp = [H_temp,ortho_func(coll_pts(i,j))];
        end
    end
    K(i,:) = H_temp;
    fprintf('Create CM equations... %d%%\n',round(i/size(coll_pts,1)*100));
end


% ------------------
subplot(2,1,1)
plot(y_fem,'b-');
subplot(2,1,2)
hold on
plot(coll_pts(:,1),'r--');
plot(coll_pts(:,2),'g--');
plot(coll_pts(:,3),'b--');
plot(coll_pts(:,4),'k--');
legend('ce','ck1','ck2','phie');
set(gca,'yscale','log');