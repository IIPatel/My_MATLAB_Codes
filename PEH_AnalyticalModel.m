clc;
clear all;
close all;
disp('PiezoSystems,Inc.:T226-A4-503%');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piezo material properties
Eps_0 = 8.854e-12;
Eps_T_33 = 3400*8.854e-12;
rho_p = 7500;
d31 = -186e-12;
sE11=16.5e-12;
sE12=-4.78e-12;
cE11 = 56e9;
% Calculate material properties for current configuration

%%%%%%% Plate configuration
% cE11 = sE11/(sE11^2-sE12^2);
% e31 = d31/(sE11+sE12);
%Eps_S_33 = EpsT_33-2*d31^2/(sE11+sE12);
% KS33=Eps_S_33/Eps_0

% Beam Configuration%%%%%%%%
% cE11 = 1/sE11;
% e31 = d31/sE11;
%e31 = -14; % fit to frequency spacing
% Eps_S_33 = Eps_T_33-d31^2/sE11;
%Eps_S_33= Eps_T_33-d31*e31;
%KS33 = Eps_S_33/Eps_0;

% Parameter assignment
c_E = cE11;
e_m = e31;
Eps_S = Eps_S_33;

% Material properties for layers
E_pzt_p = c_E; % PZT plate stiffness
rho_pzt = rho_p;

% structure
E_s = 100e9; %calculated
rho_s = 7165; %measured
nu_s = ; % beam configuration
E_s_p = E_s/(1-nu_s^2); % stiffness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define device parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 55e-3; % length in m
b = 31.8e-3; % beam width in m
%Specify layer thicknesses
t_pzt = 270e-6;
t_s = 140e-6
% Define parameters for mass at tip, not applicable for this case
L0 = 0; % length of the mass == 0 (no mass)
H0 = 70e-6; %thickness of the mass
b0 = b; % width of the mass
rho_0 = 1200; % density of mass

% mechanical damping
Zeta_m = 0.0178; % measured
ddW_b = 2.5; % input vibration amplitude
% index of pzt element above the neutral axis
ind_pzt = 3;
% Specify geometry
t_1 = [t_pzt t_s t_pzt];
E_1 = [E_pzt_p E_s_p E_pzt_p];
rho_l= [rho_pzt rho_s rho_pzt];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate neutral axis, I, EI, rho, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the neutral axis (NA)
for ii = 1: length(t_l)
       if ii == 1
            zb_l(ii) = tl(ii)/2;
       else
           zb_l(ii) = zb_l(ii-1) + t_l(ii-1)/2 + t_l(ii)/2;
end % if
EAzb(ii) = zb_l(ii)*tjl(ii)*b*EJl(ii);
EA(ii) = t_l(ii)*b*EJl(ii);
end % for
zb = sum(EAzb)/sum(EA);

% useful values for later calculation
zp_max = zb_l(ind_pzt)+t_pzt/2-zb; % max piezo distance from NA
zp_ave = abs(zb_l(ind_pzt)-zb); % distance: NA to piezo center
t_t = sum(t_l); % total modeled thickness
%determine the moment of inertia for each layer
for ii = 1:length(t_l)
I_l(ii) = 1/12*b*t_l(ii)^3 + t_l(ii)*b*(zb_l(ii)-zb)^2;
end % for ii
%Moment of area
I_s = sum(I_l);
% effective EI
EI_s = 0;
for ii=1:length(t_l)
EI_s = EI_s + E_l(ii)*I_l(ii);
end % for ii
%effective young's modulus
E_eff = EI_s/I_s;
%effective density
rho_c = 0;
for ii = 1:length(t-l)
rho_c = rho_c +tl(ii)*rhojl(ii);
end % for ii
rho_c = rho_c/t_t;
%mass per length
m = rho_c*b*t_t; % mass per length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate parameters for modal analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effect of non-point loaded proof mass
m0 = b0*(rho_c*t_t + rho_0*H0);
CG_m =(rho_c*t_t*b0*zb+rho_0*H0*b0*(t_t+H0/2))/(rho_c*t_t+rho_0*H0);
ox = L0/2; % depends on the beam config
oy = CG_m-zb; % distance from CG of mass to neutral axis

M0 = m0*L0;
M0_ = M0/(m*L);
% determine the static moment
S0 = M0*ox;
S0_ = S0/(m*L^2);
% determine the moment of inertia of the mass
for ii=1:length(t_l)
Iz(ii) = rho_l(ii)*b0*(1/12*L0*t_l(ii)^3+t_l(ii)*L0*(zb_l(ii)-CG_m)^2);
Ix(ii) = rho_l(ii)*b0*(1/12*t_l(ii)*L0^3);
end IzM = rho_0*b0*(1/12*L0*H0^3 + H0*L0*(t_t+H0/2-CG_m)^2);
IxM = rho_0*b0*(1/12*H0*L0^3);
Iyy = sum(Iz) + sum(Ix) + IzM + IxM;
I0 = Iyy+M0*(ox^2+oy^2);
I0_ = I0/(m*L^3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modal Analysis (for beam with non-point loaded proof mass)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the value for sigmaN*L - range and search increments
x = [1.87:0.0001:1.88]; % w/o mass
for ii=1:length(x)
A11(ii) = (sinh(x(ii))+sin(x(ii)))+I0_*x(ii)^3*(-cosh(x(ii))+cos(x(ii))) + S0_*x(ii)^2*(-sinh(x(ii))+sin(x(ii)));
A12(ii) = (cosh(x(ii))+cos(x(ii)))+I0_*x(ii)^3*(-sinh(x(ii))-sin(x(ii))) + S0_*x(ii)^2*(-cosh(x(ii))+cos(x(ii)));
A21(ii) = (cosh(x(ii))+cos(x(ii)))+M0_*x(ii)^l*( sinh(x(ii))-sin(x(ii))) + S0_*x(ii)^2*( cosh(x(ii))-cos(x(ii)));
A22(ii) = (sinh(x(ii))-sin(x(ii)))+M0_*x(ii)-1*( cosh(x(ii))-cos(x(ii))) + S0_*x(ii)^2*( sinh(x(ii))+sin(x(ii)));
% Calculate determinant of matrix
DET(ii) = A11(ii)*A22(ii)-A12(ii)*A21(ii);
end
% Find roots
if sign(DET(1)) == 1
Ind = find(DET<O);
elseif sign(DET(1)) == -1
    Ind = find(DET>O);
end
Ind = Ind(1);
B1_L = x(Ind); % lambda_N*L
clear x;
% calculate model parameters
B1 = B1_L/L; % lambda_N for current configuration
c1 = 1; % arbitrary constant - cancels
s1 = A12(Ind)/A11(Ind);
% natural freq from modal analysis
om_n_a = B1^2*sqrt(EI_s/m);
f  = om_n_a/2/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate lumped element model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving intergrals analytically with symbolic toolbox, following
%analytical expressions are obtained and solved for each case. Note
% these functions have been truncated!
xa = 0;
xb = L;
Int_Phi_2 = -1/8*c1^2*(-exp(-2*xb*B1)*s1^2*exp(4*xb*B1)+(truncated)
Int_ddPhi_2 = 1/8*c1^2*B1^3*exp(-2*xa*B1)+1/2*c1^2*B1^3*(truncated)
Int_Phi = c1*(sinh(xb*B1)-sin(xb*B1)-s1*cosh(xb*B1)-(truncated)
Int_ddPhi= c1*B1*sinh(xb*B1)+cl*B1*sin(xb*B1)-c1*B1*s1*(truncated)
clear xa xb;
% Mode shape value calculations
Phi_L = c1*((cosh(B1*L)-cos(B1*L))-s1*(sinh(B1*L)-sin(B1*L)));
dPhi_L = c1*B1*((sinh(B1*L)+sin(B1*L))-s1*(cosh(B1*L)-cos(B1*L)));
ddPhi_0 =c1*B1^2*((cosh(B1*0)+cos(B1*0))-s1*(sinh(B1*0)+sin(B1*0)));

%calculate coefficient matrices
% Mass matrix
for ii = 1:length(t_l)
M_1(ii) = b*t_1(ii)*rho_l(ii)*Int_Phi_2;
end % for ii
M1 = sum(M_1) + M0*Phi_L^2 + 2*S0*Phi_L*dPhi_L + I0*dPhi_L^2;

% Stiffness
for ii = 1:length(t_l)
K_l(ii) = E_l(ii)*I_l(ii)*Int_ddPhi_2;
end % for ii
K1 = sum(K_l);

% Coupling matrix
Theta = b*zp_ave*e_m*dPhi_L; % for a single element
Theta = Theta; % effective coupling for series connection
% Capacitance
Cp = b*L*Eps_S/t_pzt; % For a single element
Cp = Cp/2; % effective capacitance for series connection
% System coupling, kappa
Ke_2 = Theta^2/(K1*Cp);
% Input matrix
Bf = m*Int_Tsi + M0*Phi_L + 1/2*M0*dPhi_L*(2*L+L0);
% estimate the natural frequency acc to lumped model
om_n = sqrt(K1/M1); % check to compare against modal analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify frequencies and loads at which system response is analyzed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
om = [70*2*pi 95*2*pi om_n om_n*sqrt(1+Ke_2) 125*2*pi 150*2*pi]';
OM = om/om_n; % normalize
R1= [4.61 11.91 19.90 55.9 100.1 156.2 200.4].*1000;
Re = R1*Cp*om_n; % dimensionless time constant, alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SCalculate system response at these loads and frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate optimal resistances at OC and SC
OM_oc = sqrt(1+Ke_2); % Anti-resonance frequency
r_opt = (GM_oc^4 + (4*Zeta_m^2-2)*OM_oc^2 + 1)/(OM_oc^6 +(4*Zeta_m^2-2*(1+Ke_2))*OM_oc^4+(1+Ke_2)^2*OM_oc^2);
Re_opt_oc = sqrt(r_opt);
R1_opt_oc = round(Re_opt_oc/(Cp*om_n));

OM_sc = 1; % Resonance frequency
r_opt = (OM_sc^4 + (4*Zeta_m^2-2)*OM_sc^2 + 1)/(OM_sc^6 + (4*Zeta_m^2-2*(1+Ke_2))*OM_sc^4 + (1+Ke_2)^2*OM_sc^2);
Re_opt_sc = sqrt(r_opt);
R1_opt_sc = round(Re_opt_sc/(Cp*om_n));
R1_opt_df = R1_opt_oc-R1_opt_sc;

% Calculate the optimum frequency ratio for a given electrical load
factor_H = 4*Zeta_m^2 - 2*(1+Theta^2/(K1*Cp));
% calculate optimum frequency numerically
for ii=1:length(Re)
    OM_opt_roots = roots([2*Re(ii)^2 (1+factor_H*Re(ii)^2) 0 -1]);
    for jj=1:length(OM_opt_roots)
        if real(OM_opt_roots(jj))>0
                OM_opt(ii) = sqrt(OM_opt_roots(jj));
        end
    end
end
% Use coefficients to calculate the response of the system
for ii=1:length(Re)
    for jj=1:length(OM)
        den = sqrt((1-(1+2*Zeta_m*Re(ii))*OM(jj)^2)^2+((2*Zeta_m+(1+Theta^2/(K1*Cp))*Re(ii))*OM(jj)-Re(ii)*OM(jj)^3)^2);
        Z(ii,jj) = (Bf/K1*sqrt(1+Re(ii)^2*OM(jj)^2))*ddW_b/den;
        W_tip(ii,jj) = Phi_L*Z(ii,jj);
        Rot_tip(ii,jj) = dPhi_L*Z(ii,jj);
        Str_m(ii,jj)= -zp_max*ddPhi_0*Z(ii,jj);
        Vp(ii,jj) = Bf*Re(ii)*Ke_2*OM(jj))*ddW_b/abs(Theta)/den;
        Pout(ii,jj) = om_n/K1*Re(ii)*Ke_2*OM(jj)^2 / den^2;
        Ip(ii,jj) = Pout(ii,jj)/Vp(ii,jj);
    end
end
% Calculate response for power-optimal resistance at each frequency
for ii=1:length(OM)
    Re_opt_f (ii)=sqrt((OM(ii)^4+(4*Zeta_m^2-2)*OM(ii)^2+1)/(OM(ii)^6+(4*Zeta_m^2-2*(+Ke_2))*M(ii)^4+(1+Ke_2)^2*OM(ii)^2));
    R1_opt_f(ii) = Re_opt_f(ii)/(Cp*om_n);
    den_opt=sqrt((1-(    1+2*Zeta_m*Re_opt_f(ii))*OM(ii)^2)^2+((2*Zeta_m+(1+Theta^2/(K1*Cp))*Re_opt_f (ii))*OM(ii)-Re_opt_f(ii)*OM(ii)^3)^2);
    Z_opt(ii)= Bf/K1*sqrt(1+(Re_opt_f(ii)*OM(ii))^2)*ddW_b/den_opt;
    W_tip_opt(ii) = Phi_L*Z_opt(ii);
    Rot_tip_opt(ii) = dPhi_L*Z_opt(ii);
    Str_m_opt(ii) = -zpmax*ddPhi_0*Z_opt(ii);
    Vp_opt(ii) = Bf*Re_opt_f(ii)*Ke_2*M(jj))*ddW_b/abs(Theta)/den;
    Pout_opt(ii)= om_n/K*Re_opt_f(ii)*Ke_2*OM(jj)^2 / den^2;
    I_opt(ii) = Pout_opt(ii)/Vp_opt(ii);
end
% calculate mechanical response across length of structure
x_a = [0:L/100:L]';
for ii=1:length(x_a)
    den = sqrt((1-(1+2*Zeta_m*Re_opt_sc)*OM_sc^2)^2+((2*Zeta_m+(1+Ke_2)*Re_opt_sc)*OM_sc-Re_opt_sc*OM_sc^3)^2);
    Z_1(ii) = Bf/K1*sqrt(1+(Re_opt_sc*OM_sc)^2)*ddW_b/den
    W_1(ii) = cl*((cosh(B1*x_a(ii))-cos(B1*x_a(ii)))-s1*(sinh(B1*x_a(ii))-sin(B1*x_a(ii))))*Z_1(ii);
    Str_1(ii) = -zp-max*(c1*B^2*((cosh(B*x_a(ii))+cos(B1*x_a(ii)))-s1*(sinh(B1*x_a(ii))+sin(B1*xv(ii)))))*Z_1(ii);
end
% calculate device power density
% operating volume
volume = (2*max(abs(W_tip_opt))+H0*cos(max(abs(Rot_tip_opt)))+L0*sin(max(abs(Rot_tip_opt))))*(L+L0)*(max(b,b0))*1e6;
P_dens = (max(Pout_opt(:))/1e-6)/volume;

