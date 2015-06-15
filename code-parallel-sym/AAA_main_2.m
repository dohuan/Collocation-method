% AAA_main.m : the main file for G&R of an AAA
% 4 collagen fiber families, SM, and elastin 
% solid output file   name+"_s%i.dat" %i=0, 1, 2, 3, ....
% for testing code use >> AAA_main(0.01, 0.01, 1000, 'test')
% Latest update (March 17, 2010) 


function [] = AAA_main(k_sigma_f, k_sigma_m, days, name, parallel, damage_params)

% global kc P_a r_h H_h nu_e0 nu_f0 nu_m0 phi0 G_h G_e G_m Sa La_M La_0 sigma_f0 sigma_m0 n_elm kq_c kq_m age_max op_time

%addpath('Utility');

if nargin < 5
    parallel = false;
end

if nargin < 6
    damage_params.mu = 15;
    damage_params.sigma = 5;
    damage_params.k = .5;
end

format long e;

n_elm=300;             %number of element
n_dt = 1/5;             %time steps in a day
Length = 30;           %longitudinal length =Length *r_h; 
age_max = 300;  % maximum life span of collagen and smooth muscle (days)
kq_c=1/50.0;  % the rate of degradation for collagen
kq_m=1/50.0;  %  .. for smooth muscle

op_time = 210000;

%parameters for constitutive functions (these should be obtained from
%parameter estimation with a healthy aorta
%also note that these parameter with the deposition stretch should satisfy
%the homeostatic stress
%kc=[c1 (elastin), c2(fiber), c3(fiber), c4 (SM), c5 (SM)]
% axis 1 -circumferential direction
% axis 2 -longitudinal direction

rho=1050;  %density of the wall
kc = [1.4704e+002  3.1804e+003  2.1766e+001  8.3965e+001  5.4094e+000]; 

P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
nu_e0 = 0.2;                       %mass fraction of elastin at an initial state
nu_m0 = 0.2;                      %mass fraction of SM
nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
                                   % nu_e0+nu_m0+sum(nu_f0)=1.0
phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)

G_h = 1.04 ;                     % homeostatic stretch of fibers
G_m = 1.1 ;
G_e = [1.086  1.0906];          % initial stretch of elastin layer

sigma_f0 = 170.35*10^3 ;               % homeostatic stress of collagen (pa)
sigma_m0 = 48.438*10^3 ;               % homeostatic stress for the sum of passive + active

%parameters for active tone
Sa = 36e+003;  % Max. vasoactive parameter (Pa)
La_M = 1.4 ;
La_0 = 0.8 ;

% checking the compatibility conditions between parameters, the
% homeostatic stresses, and balance equation
t_11_c = rho*nu_f0(2)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h ...
    + 2* rho*nu_f0(3)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h*sin(phi0)^2;
t_22_c = rho*nu_f0(1)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h ...
    + 2* rho*nu_f0(3)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h*cos(phi0)^2;
t_m = rho*dWmdx(G_m, kc)*G_m+Sa*(1-(La_M-1.0)^2/(La_M-La_0)^2);
t_e =  rho*dWedx(1, G_e(1), G_e(2), kc)*G_e(1);

% Calculation of in vivo thickness

stress = nu_e0*t_e + nu_m0*t_m + (1-nu_e0-nu_m0)*t_11_c;
H_h =  1.0556e-003;

% invivo thickness (m)
% disp((t_11_c-sigma_f0)/sigma_f0);
% disp((t_m-sigma_m0)/sigma_m0);
% disp((H_h-P_a*r_h/stress)/(P_a*r_h/stress));

sigma_f0=t_11_c;
sigma_m0=t_m;
H_h = P_a*r_h/stress;


% Then, initial value at the (invivo) reference configuration are calculated
% at the Gauss points.

%fa_init(name,Data_t0); %without initial damage

if parallel
    growth_par_elmloop(damage_params, days, k_sigma_f, k_sigma_m, name, Length, n_dt, floor(100*n_dt), kc, P_a, r_h, H_h, nu_e0, nu_f0, nu_m0, phi0, G_h, G_e, G_m, Sa, La_M, La_0, sigma_f0, sigma_m0, n_elm, kq_c, kq_m, age_max, op_time);
else
    % global kc P_a r_h H_h nu_e0 nu_f0 nu_m0 phi0 G_h G_e G_m Sa La_M La_0 sigma_f0 sigma_m0 n_elm kq_c kq_m age_max op_time; %#ok<TLEV>
    growth(days, k_sigma_f, k_sigma_m, name, Length, n_dt, floor(100*n_dt));
end

% clear *

end



