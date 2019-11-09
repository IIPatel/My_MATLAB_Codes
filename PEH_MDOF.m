clc
clear all
close all

%Adapted from Mracek and du Toit's Theses for MS Aeronautics and Astronautics

%input parameters
Zeta_m = 0.005;
freqin = 1000; %input vibration frequency in Hz
om_b_ref = freqin*2*pi; % input vibration frequency in rad/s
ddW_b_ref = 9.81*sqrt(freqin*1.579e-4); %input vibration_amplitude from_ psd (valid 300<freqin<1000 Hz)

% Piezo material properties
d31= -186e-12;
d33=6.70e-010;
sE33=2.07e-011;
sE11=1.65e-011;
sE12=-4.78e-012;
sE13=-8.45e-012;
sE31=-8.45e-012;
%c_E11=sE_11/(sE_11^2-sE_12^2);

Eps_T_33=8.854187817e-12*3400;
Eps_0 = 8.854187817e-12;
% assign re1evant material properties for current configuration
% plate configuration
c_E11 = sE11/(sE11^2-sE12^2);
e31 = d31/(sE11+sE12);
Eps_S_33 = Eps_T_33-2*d31^2/(sE11+sE12);
KS33 = Eps_S_33/Eps_0;

% parameter assignment
c_E = c_E11;
%c_E = 66e9
e_m = e31;
Eps_S = Eps_S_33;

% Material properties for layers
% PZT
E_pzt_p = c_E; % PZT plate stiffness
rho_pzt = 7500;

% Copper Shim Layer
E_s = 11.2e10%Young Modulus
rho_s = 8780; %Mass Density
nu_s = 0.35;% Poisson Ratio
E_s_p = E_s/(1-nu_s^2); % plate modulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters
% Lmin = 0.1e-3;
% Lmax = 1.5e-3;
% Lrange = [Lmin:(Lmax-Lmin)/500:Lmax];
% for ll=1:length(Lrange)

L = 76.68e-3 % length in m

%L = Lrange(ll)

b = 33e-3; % beam width in m

% specify thicknesses
t_pzt = 0.22e-3;
t_s = 0.22e-3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters for mass at_the end
% L0min = 0.05e-3;
% L0max = 0.5e-3;
% L0range = [L0min: (L0max-L0min) /500:L0max];
% for lo = 1:length(L0range)
  L0 = 10e-2 % length of proof mass
  b0 = 2e-2;
  H0 = 2e-2; % height of the_mass w/o midd1e pt_layer, constrained by fabrication

rho_bf = 7850;

%Calculate taper s1ope:
taper = 0;%0.5*(b0-b)/L;
fparea = 2*0.5*L*b + b0*L; %foot print area of the beam, not inc1uding the proof mass

%Create width matrix:
xstep = L/100;
xx = (0:xstep:L);
for jj = 1:length(xx)
    bx(jj) = b + 2*taper*xx(jj);
end

%plot (xx,bx/2)
%hold on
%plot(xx,-bx/2)
%hold off
%effective fractional area of e1ectrode pads, assuming 10um each sidemof beam:

% Specify layers
t_l = [t_pzt t_s t_pzt];
E_l = [E_pzt_p E_s_p E_pzt_p];
rho_l = [rho_pzt rho_s rho_pzt];

totalthick = sum(t_l);

% determine the neutral axis
for jj = 1: length(bx) %sum over beam length
    for ii = 1: length(t_l) %sum over beam thickness
        if ii == 1
            zb_l(ii) = t_l(ii)/2;
        else
            zb_l(ii) = zb_l(ii-1) + t_l(ii-1)/2 + t_l(ii)/2;

            % zb of individual layers from_ ground
        end % if
        EAzb(ii,jj) = zb_l(ii)*t_l(ii)*bx(jj)*E_l(ii);
        EA(ii,jj) = t_l(ii)*bx(jj)*E_l(ii);
    end % for
    EAzbsum = sum(EAzb);
    EAsum = sum(EA);
    zb(jj)=EAzbsum(jj)/EAsum(jj); %note: zb is constant with x b/c the b(x) in EA divides out
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% useful parameters
ind_pzt = 1; %update if changing layers under pzt
%zp max = zb_l(ind_pzt)+t_pzt/2-zb; from_ uniform case
zp_max = max(zb_l(ind_pzt)+t_pzt/2-zb); % max piezo distance from_neutral axis set manually
%zp_ave = abs(zb_l(ind_pzt)-zb); from_ uniform case
zp_ave = abs(zb_l(ind_pzt)-mean(zb)); % distance between neutral axis and midd1e of piezo e1ement
t_t = sum(t_l);%Total modeled thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine the_mom_ent of inertia and effective EI for each layer
for jj = 1: length(bx)
    for ii = 1:length(t_l)
        I_l(ii,jj) = 1/12*bx(jj)*t_l(ii)^3 + t_l (ii) *bx(jj) * (zb_l(ii)-zb (jj))^2;
        I_sx = sum(I_l,1);
        EI_l(ii,jj) = E_l(ii)*I_l(ii,jj);
    end
    EI_sx = sum(EI_l,1);
    E_effx(jj) = EI_sx(jj)/I_sx(jj); % effective young's modulus as a function of x
end

% effective density, NOT a function of x
rho_c = 0;
for ii = 1:length(t_l)
    rho_c = rho_c+t_l(ii)*rho_l(ii);
end % for ii
rho_c = rho_c/t_t;

% mass per length as a function of x
for jj = 1:length(bx)
    m(jj) = rho_c*bx(jj)*t_t;
end

% determine the_mass per length for the Proof Mass
%proof mass assumed uniform in x
m0= b0*(rho_c*t_t + rho_bf*H0);
CG_m = 0.33e-3%(rho_c*t_t*b0*zb(length(bx)) + rho_bf*H0*b0*(t_t+H0/2))/(rho_c*t_t + rho_bf*H0);
ox = L0/2; % depends on the beam config
oy = CG_m-zb(length(bx)); % distance from_ CG of mass neutral axis

% determine proof mass effective parameters for modal analysis
M0 = 0.246615023306799;   %0.123307511653399%m0*L0;
Mbi = 0.011580133767840;
m = Mbi/L;
%M0_ = M0/(m*L);
M0_= M0/(mean(m))*L
%(from_ uniform case)
% determine the static moment
S0= M0*ox;
%S0_= S0/(m*L^2); (from_ uniform case)
S0_= S0/(mean(m)*L^2);
% determine the moment of inertia of the mass
for ii=1:length(t_l)
    Iz(ii) = rho_l(ii)*b0*(1/12*L0*t_l(ii)^3 + t_l(ii)*L0*(zb_l(ii)-CG_m)^2);
    Ix (ii) = rho_l(ii)*b0*(1/12*t_l(ii)*L0^3);
end
IzM = rho_bf*b0*(1/12*L0*H0^3 + H0*L0*(t_t+H0/2-CG_m)^2);
IxM = rho_bf*b0*(1/12*H0*L0^3);
Iyy = sum(Iz) + sum(Ix) + IzM + IxM;
I0 = 0.001444753011539;  % 0.000720064489926%Iyy+M0*(ox^2+oy^2);
%IO = IO/(m*L^3);(from_ uniform case)
I0_= I0/(mean(m))*L^3

%calculate static def1ection of beam:
for jj=1:length(xx)
    %from_ point load:
    stdefl(jj) = (M0*9.81*xx(jj)^2)/(6*EI_sx(jj))*(3*L-xx(jj));
    %from_ moment:
    stdef2(jj) = (S0*xx(jj)^2)/(2*EI_sx(jj));
    %superposition:
    stdefT (jj) =stdefl (jj) +stdef2 (jj);
end
staticdef = stdefT(length(xx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the specific configuration - modal parameters
% find the value for Beta1
%x = [0.25:0.005:0.5];
x = [0.5:0.000005:0.508]; % x = lambda bar; narrow range after successive runs

for ii=1:length(x)
    A11(ii) = (sinh(x(ii))+sin(x(ii))) + I0_*x(ii)^3*(-cosh(x(ii))+cos(x(ii))) + S0_*x(ii)^2*(-sinh(x(ii))+sin(x(ii)));
    A12(ii) = (cosh(x(ii))+cos(x(ii))) + I0_*x(ii)^3*(-sinh(x(ii))-sin(x(ii))) + S0_*x(ii)^2*(-cosh(x(ii))+cos(x(ii)));

    A21(ii) = (cosh(x(ii))+cos(x(ii))) + M0_*x(ii)^1*( sinh(x(ii))-sin(x(ii))) + S0_*x(ii)^2*( cosh(x(ii))-cos(x(ii)));
    A22(ii) = (sinh(x(ii))-sin(x(ii))) + M0_*x(ii)^1*( cosh(x(ii))-cos(x(ii))) + S0_*x(ii)^2*( sinh(x(ii))+sin(x(ii)));

    % determine the determinant
    DET(ii) = A11(ii)*A22(ii)-A12(ii)*A21(ii);
end
    if sign(DET(1)) == 1
        Ind = find(DET<0);
    elseif sign(DET(1)) == -1
        Ind = find(DET>0);
    end
    
    if numel(Ind)== 0
        pickmetric = 0;
        disp('error with Ind')
        return
    end
Inda = Ind(1);
Gamma_ = x(Inda)  % this is Beta*L
%disp(['BetaL = ', num2str(Gamma_)]);
%c1ear x;

B1 = Gamma_/L;
c1 = 1; % arbitrary constant - Cancels out when determining the displacement and other parameters
s1 = A12(Inda)/A11(Inda);
% natural freq acc to modal analysis
om_n_a = B1^2*sqrt(EI_sx/(rho_c*t_t*fparea));
f = om_n_a/2/pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate coefficient matrices:
Phi_L = c1*((cosh(B1*L)-cos(B1*L))-s1*(sinh(B1*L)-sin(B1*L)));
dPhi_L = c1*B1*((sinh(B1*L)+sin(B1*L))-s1*(cosh(B1*L)-cos(B1*L)));
ddPhi_0 =c1*B1^2*((cosh(B1*0)+cos(B1*0))-s1*(sinh(B1*0)+sin(B1*0)));

for jj = 1:length(xx)
    Tsilx(jj) = c1*( (cosh(B1*xx(jj))-cos(B1*xx(jj))) - s1*(sinh(B1*xx(jj))-sin(B1*xx(jj))));
    ddTsilx(jj) = c1*( (B1^2*cosh(B1*xx(jj))+B1^2*cos(B1*xx(jj))) - s1*(B1^2*sinh(B1*xx(jj))+B1^2*sin(B1*xx(jj))));
end

%Mass matrix:
for ii = 1:length(t_l)
    for jj = 1:length(bx)
        Mint(ii,jj) = bx(jj)*t_l(ii)*rho_l(ii)*Tsilx(jj)^2*xstep;
    end
    M_l(ii) = sum(Mint(ii,:));
end

for jj=1:length(xx)
Msum(jj) = sum(Mint(:,jj)/xstep);
end
Ml = sum(M_l) + M0*Phi_L^2 + 2*S0*Phi_L*dPhi_L + I0*dPhi_L^2; %effective_mass

%Stiffness matrix:
for jj = 1:length(xx)
    for ii = 1:length(t_l)
        Kint(ii,jj) = E_l(ii)*I_l(ii,jj)*ddTsilx(jj)^2*xstep;
    end
    Kx(jj) = sum(Kint(:,jj)/xstep);
end
K1 = sum(Kx.*xstep);

%Coupling matrix
for jj=1:length(bx)
    Thetax(jj) = bx(jj)*zp_ave*e_m*ddTsilx(jj)*xstep;
    Thetax2(jj) = Thetax(jj)/xstep;
end %for jj (length)
Theta = sum(Thetax);

%Capacitance_matrix
for jj = 1:length(bx)
    Cpx(jj) = (bx(jj)*Eps_S/t_pzt)*xstep;
    Cpx2(jj) = Cpx(jj)/xstep;
end %for jj (length)
Cp = sum(Cpx);

Ke_2 = Theta^2/(K1*Cp);

%Input matrix adjusted to account for proof mass
for jj = 1:length(xx)
    Bf_int(jj) = rho_c*t_t*bx(jj)*Tsilx(jj)*xstep;
    Bfx(jj) = Bf_int(jj)/xstep;
end
Bf = sum(Bf_int) + M0*Phi_L + 1/2*M0*dPhi_L*(2*L+L0);

% estimate the natural frequency acc to lumped model
om_n = sqrt(K1/Ml);
%om_ = om_n*sqrt(1+Ke_2);
%OM_ = om_/om_n;

OM_ = om_b_ref/om_n; %want to make sure that the OM ref1ects off optimum operation
ddW_b = ddW_b_ref;

% Analyse the system%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify frequencies and loads at which system response is analyzed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
om_=2*pi.*[0:0.001:10]';%om_ = 2*pi.*[0 40 50 70 95 125 150]';%om_n/2/pi om_n*sqrt(1+Ke_2)/2/pi
OM_ = om_/om_n; % normalize
R1= [4.61 11.91 19.90 55.9 100.1 156.2 200.4].*1000;
Re = R1*Cp*om_n; % dimensionless time constant, alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SCalculate system response at these loads and frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate optimal resistances at OC and SC
OM_oc = sqrt(1+Ke_2); % Anti-resonance frequency
r_opt = (OM_oc^4 + (4*Zeta_m^2-2)*OM_oc^2 + 1)/(OM_oc^6 +(4*Zeta_m^2-2*(1+Ke_2))*OM_oc^4+(1+Ke_2)^2*OM_oc^2);
Re_opt_oc = sqrt(r_opt);
R1_opt_oc = round(Re_opt_oc/(Cp*om_n))

OM_sc = 1; % Resonance frequency
r_opt = (OM_sc^4 + (4*Zeta_m^2-2)*OM_sc^2 + 1)/(OM_sc^6 + (4*Zeta_m^2-2*(1+Ke_2))*OM_sc^4 + (1+Ke_2)^2*OM_sc^2);
Re_opt_sc = sqrt(r_opt);
R1_opt_sc = round(Re_opt_sc/(Cp*om_n))
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
    for jj=1:length(OM_)
        den = sqrt((1-(1+2*Zeta_m*Re(ii))*OM_(jj)^2)^2+((2*Zeta_m+(1+Theta^2/(K1*Cp))*Re(ii))*OM_(jj)-Re(ii)*OM_(jj)^3)^2);
        Z(ii,jj) = (Bf/K1*sqrt(1+Re(ii)^2*OM_(jj)^2))*ddW_b/den;
        W_tip(ii,jj) = Phi_L*Z(ii,jj);
        Rot_tip(ii,jj) = dPhi_L*Z(ii,jj);
        Str_m(ii,jj)= -zp_max*ddPhi_0*Z(ii,jj);
        Vp(ii,jj) = Bf*Re(ii)*Ke_2*OM_(jj)*ddW_b/abs(Theta)/den;
        Pout(ii,jj) = om_n/K1*Re(ii)*Ke_2*OM_(jj)^2 / den^2;
        Ip(ii,jj) = Pout(ii,jj)/Vp(ii,jj);
    end
end
% Calculate response for power-optimal resistance at each frequency
for ii=1:length(OM_)
    Re_opt_f (ii)=sqrt((OM_(ii)^4+(4*Zeta_m^2-2)*OM_(ii)^2+1)/(OM_(ii)^6+(4*Zeta_m^2-2*(+Ke_2))*OM_(ii)^4+(1+Ke_2)^2*OM_(ii)^2));
    R1_opt_f(ii) = Re_opt_f(ii)/(Cp*om_n);
    den_opt=sqrt((1-(    1+2*Zeta_m*Re_opt_f(ii))*OM_(ii)^2)^2+((2*Zeta_m+(1+Theta^2/(K1*Cp))*Re_opt_f (ii))*OM_(ii)-Re_opt_f(ii)*OM_(ii)^3)^2);
    Z_opt(ii)= Bf/K1*sqrt(1+(Re_opt_f(ii)*OM_(ii))^2)*ddW_b/den_opt;
    W_tip_opt(ii) = Phi_L*Z_opt(ii);
    Rot_tip_opt(ii) = dPhi_L*Z_opt(ii);
    Str_m_opt(ii) = -zp_max*ddPhi_0*Z_opt(ii);
    Vp_opt(ii) = Bf*Re_opt_f(ii)*Ke_2*OM_(jj)*ddW_b/abs(Theta)/den;
    Pout_opt(ii)= om_n/K1*Re_opt_f(ii)*Ke_2*OM_(jj)^2 / den^2;
    I_opt(ii) = Pout_opt(ii)/Vp_opt(ii);
end
% calculate mechanical response across length of structure
x_a = [0:L/100:L]';
for ii=1:length(x_a)
    den = sqrt((1-(1+2*Zeta_m*Re_opt_sc)*OM_sc^2)^2+((2*Zeta_m+(1+Ke_2)*Re_opt_sc)*OM_sc-Re_opt_sc*OM_sc^3)^2);
    Z_1(ii) = Bf/K1*sqrt(1+(Re_opt_sc*OM_sc)^2)*ddW_b/den;
    W_1(ii) = c1*((cosh(B1*x_a(ii))-cos(B1*x_a(ii)))-s1*(sinh(B1*x_a(ii))-sin(B1*x_a(ii))))*Z_1(ii);
    Str_1(ii) = -zp_max*(c1*B1^2*((cosh(B1*x_a(ii))+cos(B1*x_a(ii)))-s1*(sinh(B1*x_a(ii))+sin(B1*x_a(ii)))))*Z_1(ii);
end
% calculate device power density
% operating volume
volume = (2*max(abs(W_tip_opt))+H0*cos(max(abs(Rot_tip_opt)))+L0*sin(max(abs(Rot_tip_opt))))*(L+L0)*(max(b,b0))*1e6;
P_dens = (max(Pout_opt(:))/1e-6)/volume;


% calculate the sc andoc resonance freq ratios%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OM_oc = sqrt(1+Theta^2/(K1*Cp));
% r_opt_ = (OM_oc^4 + (4*Zeta_m^2-2)*OM_oc^2 + 1)/(OM_oc^6 +(4*Zeta_m^2-2*(1+Theta^2/(K1*Cp)))*OM_oc^4 + (1+Theta^2/(K1*Cp))^2*OM_oc^2);
% Re_opt_oc = sqrt(r_opt_);
% Rl_opt_oc = round(Re_opt_oc/(Cp*om_n));
% OM_sc = 1;
% r_opt_ = (OM_sc^4 + (4*Zeta_m^2-2)*OM_sc^2 + 1)/(OM_sc^6 +(4*Zeta_m^2-2*(1+Theta^2/(K1*Cp)))*OM_sc^4 + (1+Theta^2/(K1*Cp) )^ 2*OM_sc^2);
% Re_opt_sc = sqrt(r_opt_);
% Rl_opt_sc = round(Re_opt_sc/(Cp*om_n));
% Rl_opt_df = Rl_opt_oc-Rl_opt_sc;
% 
% % Specify open circuit optimum resistance:
% Rl = Rl_opt_oc;
% Re = Rl*Cp*om_n;
% 
% % Next calculate the response of the system
% den = sqrt( (1 - (1 + 2*Zeta_m*Re)*OM_^2 )^2 + ( (2*Zeta_m + (1+Theta^2/(K1*Cp))*Re)*OM_ - Re*OM_^3 )^2 );
% Z = (Bf/K1*sqrt(1+Re^2*OM_^2))*ddW_b / den;
% W_tip = Phi_L*Z;
% Rot_tip = dPhi_L*Z;
% Str_m = -zp_max*ddPhi_0*Z;
% Vp = (Bf/(K1*Cp)*Re*Theta*OM_)*ddW_b / den;
% Pout = ((Bf^2/K1^2)*Theta^2*OM_^2*om_n*(1/Cp)*Re)*ddW_b^2 / den^2;
% Ip = Pout/Vp;
% 
% Re_opt_f = sqrt((OM_^4 + (4*Zeta_m^2-2)*OM_^2 + 1)/(OM_^6 + (4*Zeta_m^2-2*(1+Theta^2/(K1*Cp)))*OM_^4 + (1+Theta^2/(K1*Cp) )^ 2*OM_^2));
% R1_opt_f = Re_opt_f/(Cp*om_n);
% den_opt_ = sqrt( (1 - (1 + 2*Zeta_m*Re_opt_f)*OM_^2 )^2 +((2*Zeta_m + (1+Theta^2/(K1*Cp))*Re_opt_f)*OM_ - Re_opt_f*OM_^3 )^2 );
% Z_opt_ = Bf/K1*sqrt(1+(Re_opt_f*OM_)^2)*ddW_b / den_opt_;
% W_tip_opt_ = Phi_L*Z_opt_;
% Rot_tip_opt_ = dPhi_L*Z_opt_;
% Str_m_opt_ = -zp_max*ddPhi_0*Z_opt_;
% Vp_opt_ = Bf/(K1*Cp)*Re_opt_f*Theta*OM_*ddW_b / den_opt_;
% Pout_opt_ = ((Bf/K1*Theta*OM_*ddW_b)^2*om_n*(1/Cp)*Re_opt_f) / den_opt_^2;
% I_opt_ = Pout_opt_/Vp_opt_;
% 
% %output data:
opvolume = (2*max(abs(W_tip_opt)) + H0*cos(max(abs(Rot_tip_opt)))+L0*sin(max(abs(Rot_tip_opt))))*(L+L0)*(max(b,b0))*1e6; % approximate volume in cmA3
%volume = 2*max(abs(W_tip_opt_))*L*(max(b,b0))*1e6; % approximate volume in cm^3

P_dens_op = (Pout_opt/1e-6)/opvolume;
volmm3 = opvolume*1e3;
volstatic = (staticdef + t_t + H0) * (L + L0) * max(bx) *1e6;
footprint = (L+L0)*(max(b,b0))*1e6;
P_dens_st = (max(Pout_opt(:))/1e-6)/volstatic;
beammass = rho_c*mean(bx)*t_t*L ;
pmmass = rho_bf*H0*b0*L0+ rho_c*bx(length(bx))*t_t*L0;
devicemass = rho_c*mean(bx)*t_t*L + rho_c*bx(length(bx))*t_t*L0 + rho_bf*H0*b0*L0;

% om_r(ll,lo) = om_n/2/pi;
% om_ar(ll,lo) = om_n/2/pi*OM_oc;
% OM_ocout=OM_oc;
% Pout_out(ll,lo) =(Pout_opt_)*1e6;
% Popdens_(ll,lo) =P_dens_op;
% Pstdens_(ll,lo) =P_dens_st;
% SpecificPower(ll,lo) = Pout_opt_/devicemass;
% maxustrain(ll,lo)=abs(Str_m_opt_)/1e-6;
% maxV(ll,lo)=abs(Vp_opt_);
% maxuA(ll,lo)=abs(I_opt_)/1e-6;
% wtip(ll,lo)=abs(W_tip_opt_)*1e6;

subplot(221)
plot(om_/2/pi, Pout(1,:))
xlabel('Input Frequency [Hz]')
ylabel('Output Power [W]')
subplot(222)
plot(om_/2/pi, W_tip(1,:))
xlabel('Input Frequency [Hz]')
ylabel('Tip Displacement [cm]')
subplot(223)
plot(om_/2/pi, Vp(1,:))
xlabel('Input Frequency [Hz]')
ylabel('Output Voltage across R_L =  k\Omega [mV]')
