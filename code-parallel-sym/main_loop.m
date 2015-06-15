% function [ TD_L1_part, TD_L2_part ] = main_loop( elm, k, g )
% %MAIN_LOOP Summary of this function goes here
% %   Detailed explanation goes here
% 
% % REQUIRED GLOBAL VALUES IN STRUCT
% % These do NOT change during the loop
% % n_Gauss_pt
% % Gauss_pt
% % l_elm
% % npe
% % Data
% 
% % count=count+1;
% count = ((elm-1)*g.n_Gauss_pt) + k;
% 
% Z=((elm-1)+(g.Gauss_Pt(g.n_Gauss_pt,k)+1)/2.0)*g.l_elm;
% 
% R=0.0;
% dR=0.0;
% 
% z= 0.0;
% r=0.0;
% dz =0.0;
% dr = 0.0;
% 
% for i=1:g.npe
%     index=(elm-1)*(g.npe-1)+i;
%     R = R+Data(index,2)*shp1d(index,0,Z);
%     dR = dR+Data(index,2)*shp1d(index,1,Z);
% 
%     z = z +x_np(index*2)*shp1d(index,0,Z);
%     r = r +x_np(index*2-1)*shp1d(index,0,Z);
%     dz =dz+x_np(index*2)*shp1d(index,1,Z);
%     dr = dr+x_np(index*2-1)*shp1d(index,1,Z);
% 
% end
% 
% sz2r2=sqrt(dz^2+dr^2);
% s1R2=sqrt(1+dR^2);
% 
% L1=r/R;                                                %lambda_1 (circumflex)
% L2=sz2r2/s1R2;     %lambda_theta (z-direction)
% TD_L1_part(k, nt+1) = L1;
% TD_L2_part(k, nt+1) = L2;
% 
% dL1dr=1/R;
% dL2dr2=dr/(s1R2*sz2r2);
% dL2dz2=dz/(s1R2*sz2r2);
% ddL2ddr2=dz^2 / (s1R2*sz2r2^3);
% ddL2ddz2=dr^2 / (s1R2*sz2r2^3);
% ddL2dr2dz2 = - dr*dz / (s1R2*sz2r2^3);
% 
% if nt+1 <= num_DL  %If there still remains initial mass
%     % phi = phi0
%     L_k3_parallel=sqrt( (L1*sin(phi0))^2+(L2*cos(phi0))^2);
%     L_k4_parallel=sqrt( (L1*sin(pi-phi0))^2+(L2*cos(pi-phi0))^2);
% 
%     dL_k3dL1=L1*sin(phi0)^2 / L_k3_parallel;  %dL_k4dL1=dL_k3dL1
%     dL_k3dL2=L2*cos(phi0)^2 / L_k3_parallel;
%     ddL_k3ddL1=(L2*sin(phi0)*cos(phi0))^2 /L_k3_parallel^3;
%     ddL_k3ddL2=(L1*sin(phi0)*cos(phi0))^2 /L_k3_parallel^3;
%     ddL_k3dL1dL2= - L1*L2*(sin(phi0)*cos(phi0))^2 / L_k3_parallel^3;
%     dL_k4dL1=L1*sin(pi-phi0)^2 / L_k4_parallel;
%     dL_k4dL2=L2*cos(pi-phi0)^2 / L_k4_parallel;
%     ddL_k4ddL1=(L2*sin(pi-phi0)*cos(pi-phi0))^2 /L_k4_parallel^3;
%     ddL_k4ddL2=(L1*sin(pi-phi0)*cos(pi-phi0))^2 /L_k4_parallel^3;
%     ddL_k4dL1dL2= - L1*L2*(sin(pi-phi0)*cos(pi-phi0))^2 / L_k4_parallel^3;
% 
%     Ln_k1 = G_h*L1;
%     Ln_k2 = G_h*L2;
%     Ln_k3 = G_h*L_k3_parallel;
%     Ln_k4 = G_h*L_k4_parallel;
% 
%     Mf_r=Mf0*DQ_c(nt+1)*Rm_exp(1, Z, z_0, current_t, init_dmg_t);
%     Mt = Mf_r(1)+Mf_r(2)+Mf_r(3)+Mf_r(4);
% 
%     dwdL_k1= A_11*Mf_r(1)*G_h*dWcdx(Ln_k1);
%     dwdL_k2= A_11*Mf_r(2)*G_h*dWcdx(Ln_k2);
%     dwdL_k3= A_11*Mf_r(3)*G_h*dWcdx(Ln_k3);
%     dwdL_k4= A_11*Mf_r(4)*G_h*dWcdx(Ln_k4);
% 
% 
%     ddwddL_k1=A_11*Mf_r(1)*G_h^2*ddWcddx(Ln_k1);
%     ddwddL_k2=A_11*Mf_r(2)*G_h^2*ddWcddx(Ln_k2);
%     ddwddL_k3=A_11*Mf_r(3)*G_h^2*ddWcddx(Ln_k3);
%     ddwddL_k4=A_11*Mf_r(4)*G_h^2*ddWcddx(Ln_k4);
% 
%     dwdL1=dwdL_k1+dwdL_k3*dL_k3dL1+dwdL_k4*dL_k4dL1;
%     dwdL2=dwdL_k2+dwdL_k3*dL_k3dL2+dwdL_k4*dL_k4dL2;
%     ddwddL1=ddwddL_k1+ddwddL_k3*dL_k3dL1^2+dwdL_k3*ddL_k3ddL1 + ...
%         ddwddL_k4*dL_k4dL1^2+dwdL_k4*ddL_k4ddL1;
%     ddwddL2=ddwddL_k2+ddwddL_k3*dL_k3dL2^2+dwdL_k3*ddL_k3ddL2 + ...
%         ddwddL_k4*dL_k4dL2^2+dwdL_k4*ddL_k4ddL2;
%     ddwdL1dL2=ddwddL_k3*dL_k3dL1*dL_k3dL2+dwdL_k3*ddL_k3dL1dL2+ ...
%         ddwddL_k4*dL_k4dL1*dL_k4dL2+dwdL_k4*ddL_k4dL1dL2;
% 
%     % we assume that the initial SM removed with the elastin
%     % degradation
% 
%     Mm = Mm0*DQ_m(nt+1)*Rm_exp(0, Z, z_0, current_t, init_dmg_t);
%     Mt = Mt+Mm;
%     Ln_m = G_m * L1;
%     dwdL1 = dwdL1+ A_11*Mm*G_m*dWmdx(Ln_m);
%     ddwddL1 = ddwddL1 + A_11*Mm*G_m^2*ddWmddx(Ln_m);
% 
% else  % All of initial fiber families and SM are removed
%     dwdL1=0.0;
%     dwdL2=0.0;
%     ddwddL1=0.0;
%     ddwddL2=0.0;
%     ddwdL1dL2=0.0;
%     Mt = 0.0;
%     Mm = 0.0;
% end
% 
% % strain energy due to elastin layer
% 
% Me = Me0*Rm_exp(0, Z, z_0, current_t, init_dmg_t);
% Mt = Mt+Me;
% 
% Ln_e1 = G_e(1)*L1;
% Ln_e2 = G_e(2)*L2;
% 
% dwdL_e1= A_11*Me*G_e(1)*dWedx(1,Ln_e1, Ln_e2);
% dwdL_e2= A_11*Me*G_e(2)*dWedx(2,Ln_e1, Ln_e2);
% ddwddL_e1=A_11*Me*(G_e(1))^2*ddWeddx(1,Ln_e1,Ln_e2);
% ddwddL_e2=A_11*Me*(G_e(2))^2*ddWeddx(2,Ln_e1,Ln_e2);
% ddwddL_e1dL_e2 = A_11*Me*G_e(1)*G_e(2)*ddWeddx(3,Ln_e1,Ln_e2);
% 
% dwdL1 = dwdL1 +dwdL_e1;
% dwdL2 = dwdL2 +dwdL_e2;
% ddwddL1 = ddwddL1 +ddwddL_e1;
% ddwddL2 = ddwddL2 +ddwddL_e2;
% ddwdL1dL2 = ddwdL1dL2 +ddwddL_e1dL_e2;
% 
% %-----------------------------------------------------------
% %       numerical integration
% % int(0 to t) [ ff(tau)]dtau = ff(0)*0.5*dtau+
% % ff(dtau)*dtau+ff(2*dtau)*dtau+ . .+ff(t)*0.5*dtau
% %-----------------------------------------------------------
% if nt+1 <= num_DL
%     n_tau0 = 1;
% else
%     n_tau0= nt-num_DL+2;
% end
% for n_tau=n_tau0: (nt+1)
%     tau=(n_tau-1)*dt;
%     L1_tau=TD_L1_part(k,n_tau);
%     L2_tau=TD_L2_part(k,n_tau);
% 
%     phi = phi0;  % phi(\tau) at the reference configuration
% 
%     L_k3_parallel=sqrt( (L1*sin(phi))^2+(L2*cos(phi))^2);
%     L_k4_parallel=sqrt( (L1*sin(pi-phi))^2+(L2*cos(pi-phi))^2);
% 
%     dL_k3dL1=L1*sin(phi)^2 / L_k3_parallel;  %dL_k4dL1=dL_k3dL1
%     dL_k3dL2=L2*cos(phi)^2 / L_k3_parallel;
%     ddL_k3ddL1=(L2*sin(phi)*cos(phi))^2 /L_k3_parallel^3;
%     ddL_k3ddL2=(L1*sin(phi)*cos(phi))^2 /L_k3_parallel^3;
%     ddL_k3dL1dL2= - L1*L2*(sin(phi)*cos(phi))^2 / L_k3_parallel^3;
%     dL_k4dL1=L1*sin(pi-phi)^2 / L_k4_parallel;
%     dL_k4dL2=L2*cos(pi-phi)^2 / L_k4_parallel;
%     ddL_k4ddL1=(L2*sin(pi-phi)*cos(pi-phi))^2 /L_k4_parallel^3;
%     ddL_k4ddL2=(L1*sin(pi-phi)*cos(pi-phi))^2 /L_k4_parallel^3;
%     ddL_k4dL1dL2= - L1*L2*(sin(pi-phi)*cos(pi-phi))^2 / L_k4_parallel^3;
% 
%     L_k3_tau=sqrt( (L1_tau*sin(phi))^2+(L2_tau*cos(phi))^2);
%     L_k4_tau = L_k3_tau;
% 
%     if ((n_tau==n_tau0) || (n_tau==nt+1))
%         w_dt=dt*0.5;
%     else
%         w_dt=dt;
%     end
% 
%     mc_tau=[TD_mc1(count,n_tau), TD_mc2(count, n_tau), ...
%         TD_mc3(count, n_tau), TD_mc3(count, n_tau)]*q_i(0, current_t-tau);
%     mm_tau = TD_mm(count, n_tau) * q_i(1, current_t-tau)...
%         *Rm_exp(0, Z, z_0, current_t, init_dmg_t);
%     Mm = Mm+mm_tau*w_dt;
%     Mt = Mt+ (mc_tau(1)+mc_tau(2)+mc_tau(3)+mc_tau(4))*w_dt...
%         +mm_tau*w_dt;
% 
%     Ln_k1 = G_h*L1 / L1_tau;
%     Ln_k2 = G_h*L2 / L2_tau;
%     Ln_k3 = G_h*L_k3_parallel / L_k3_tau;
%     Ln_k4 = G_h*L_k4_parallel / L_k4_tau;
%     Ln_m =  G_m*L1/L1_tau;
% 
%     dwdL_k1 = A_11* mc_tau(1) * (G_h/L1_tau)*dWcdx(Ln_k1);
%     dwdL_k2 = A_11* mc_tau(2) * (G_h/L2_tau)*dWcdx(Ln_k2);
%     dwdL_k3 = A_11* mc_tau(3) * (G_h/L_k3_tau)*dWcdx(Ln_k3);
%     dwdL_k4 = A_11* mc_tau(4) * (G_h/L_k4_tau)*dWcdx(Ln_k4);
% 
%     ddwddL_k1 = A_11*mc_tau(1) * (G_h/L1_tau)^2*ddWcddx(Ln_k1);
%     ddwddL_k2 = A_11*mc_tau(2) * (G_h/L2_tau)^2*ddWcddx(Ln_k2);
%     ddwddL_k3 = A_11*mc_tau(3) * (G_h/L_k3_tau)^2*ddWcddx(Ln_k3);
%     ddwddL_k4 = A_11*mc_tau(4) * (G_h/L_k4_tau)^2*ddWcddx(Ln_k4);
% 
%     dwdL1=dwdL1 + (dwdL_k1+dwdL_k3*dL_k3dL1+dwdL_k4*dL_k4dL1)*w_dt;
%     dwdL2=dwdL2 + (dwdL_k2+dwdL_k3*dL_k3dL2+dwdL_k4*dL_k4dL2)*w_dt;
%     ddwddL1=ddwddL1+ (ddwddL_k1+ddwddL_k3*dL_k3dL1^2+dwdL_k3*ddL_k3ddL1 + ...
%         ddwddL_k4*dL_k4dL1^2+dwdL_k4*ddL_k4ddL1)*w_dt ;
%     ddwddL2=ddwddL2+ (ddwddL_k2+ddwddL_k3*dL_k3dL2^2+dwdL_k3*ddL_k3ddL2 + ...
%         ddwddL_k4*dL_k4dL2^2+dwdL_k4*ddL_k4ddL2)*w_dt ;
%     ddwdL1dL2=ddwdL1dL2+(ddwddL_k3*dL_k3dL1*dL_k3dL2+dwdL_k3*ddL_k3dL1dL2+ ...
%         ddwddL_k4*dL_k4dL1*dL_k4dL2+dwdL_k4*ddL_k4dL1dL2) * w_dt;
% 
%     %Contribution of SM
%     dwdL1 = dwdL1+A_11* mm_tau * (G_m/L1_tau)*dWmdx(Ln_m)*w_dt;
%     ddwddL1 = ddwddL1+A_11*mm_tau * (G_m/L1_tau)^2*ddWmddx(Ln_m)*w_dt;
% end
% 
% f_act = 1-((La_M-1.0)/(La_M-La_0))^2;
% dwdL1_act = A_11*Mm*(1.0/L1)*(Sa/rho)*f_act;
% %ddwddL1_act = A_11*Mm*Sa/rho*(1/L1_s)^2*2.0*(La_M-L1/L1_s)/(La_M-La_0)^2;
% ddwddL1_act = -A_11*Mm*(1.0/L1^2)*(Sa/rho)*f_act;
% 
% dwdL1 = dwdL1+dwdL1_act;
% ddwddL1 = ddwddL1+ddwddL1_act;
% 
% dwdr =dwdL1*dL1dr;
% dwdr2=dwdL2*dL2dr2;
% dwdz2=dwdL2*dL2dz2;
% ddwddr = ddwddL1*dL1dr^2;
% ddwddr2= ddwddL2*dL2dr2^2+dwdL2* ddL2ddr2;
% ddwddz2= ddwddL2*dL2dz2^2+dwdL2*ddL2ddz2;
% ddwdrdr2=ddwdL1dL2*dL1dr*dL2dr2;
% ddwdrdz2=ddwdL1dL2*dL1dr*dL2dz2;
% ddwdr2dz2 = ddwddL2*dL2dr2*dL2dz2+dwdL2*ddL2dr2dz2;
% 
% %Global K , K_ij
% WI=Gauss_Wi(n_Gauss_pt, k)*l_elm/2;
% for i=1:npe
%     ii = (elm-1)*(npe-1)+i;
%     ii_1=ii*2-1;
%     ii_2=ii*2;
%     Si = shp1d(ii,0,Z);
%     dSi=shp1d(ii,1,Z);
% 
% 
%     glf(ii_1)=glf(ii_1)+2*pi*((dwdr*Si + dwdr2*dSi)*R*s1R2-beta*P*dz*r*Si)*WI;
%     % ii_1 = (elm-1)*(npe-1)+(1:npe) = (elm-1)
%     % add this to glf(ii_1);
% 
%     if current_t>=op_time
%         if r_stent > r
%             F_penalty = -1.0*sz2r2*k_stent*(r_stent-r)*Si;
%             glf(ii_1)=glf(ii_1)+2*pi*F_penalty*WI;
%         end
%     end
%     glf(ii_2)=glf(ii_2)+2*pi*(dwdz2*dSi*R*s1R2 + beta*P*dr*r*Si)*WI;
% 
%     for j=1:npe
%         jj= (elm-1)*(npe-1)+j;
%         jj_1=jj*2-1;
%         jj_2=jj*2;
%         SiSj=shp1d(ii,0,Z)*shp1d(jj,0,Z);
%         SidSj=shp1d(ii,0,Z)*shp1d(jj,1,Z);
%         dSiSj=shp1d(ii,1,Z)*shp1d(jj,0,Z);
%         dSidSj=shp1d(ii,1,Z)*shp1d(jj,1,Z);
% 
% 
%         glk(ii_1,jj_1)=glk(ii_1,jj_1)+2*pi*((ddwddr*SiSj+ddwdrdr2*(SidSj+dSiSj)...
%             +ddwddr2*dSidSj)*R*s1R2-beta*P*dz*SiSj)*WI;
%         glk(ii_1,jj_2)=glk(ii_1,jj_2)+2*pi*((ddwdrdz2*SidSj+ ddwdr2dz2*dSidSj)*...
%             R*s1R2-beta*P*r*SidSj)*WI;
% 
%         if current_t>=op_time
%             if r_stent > r
%                 K11_penalty = (sz2r2*k_stent*SiSj-(dr/sz2r2)*k_stent*(r_stent-r)*SidSj);
%                 K12_penalty =-1*(dz/sz2r2)*k_stent*(r_stent-r)*SidSj;
%                 glk(ii_1,jj_1)=glk(ii_1,jj_1)+2*pi*K11_penalty*WI;
%                 glk(ii_1,jj_2)=glk(ii_1,jj_2)+2*pi*K12_penalty*WI;
%             end
%         end
%         glk(ii_2,jj_2)=glk(ii_2,jj_2)+2*pi*(ddwddz2*dSidSj*R*s1R2+beta*P*r*SidSj)*WI;
% 
%         glk(ii_2,jj_1)=glk(ii_2,jj_1)+2*pi*((ddwdrdz2*dSiSj+ddwdr2dz2*dSidSj)*...
%             R*s1R2+beta*P*(dr*SiSj+r*SidSj))*WI;
% 
%     end
% end
% 
% end
% 
