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
