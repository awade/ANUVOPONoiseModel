%% NOISE BUDGET SCRIPT FOR ANU VOPO PROJECT

% Authors: Andrew Wade
% Date: 6 Nov 2015
%
% Notes: 
% -All units in standard SI unless otherwise stated
% - fluctating components writeen as x_delta
% - Steady state components written as x_ss
%
% Previous Versions: NA
%
%
% Comments: The idea is to build this script out to the point where
% indivual noise models can be spun out into indivdual functions that can
% be called from the main script.  This way new and updated models can
% easly be added/modified and swapped out.
%


clear all  %Clear the decks
close all

%Standard physical constants
c = 3e8; %[m/s]
h = 6.626e-34; %[J.s]
n = 1; %Refractive index of medium
lambdaFund = 1064e-9; %[m] wave length of light
lambdaHarm = 532e-9; %[m] wave length of light


%%%%%%%%%%%%%%%%%%%%% Cavity parameters %%%%%%%%%%%%%%%%%%%%%%%
L = 0.345+0.01*0.83; %[m] total effective cavity round trip length

%Fundamental field parameters
Rain = 0.845;%+0.0001; %Input coupler reflectivity
Raout = 1; %Output (transmission) coupler reflectivity
Lossa = 1-0.989; 0.0; %Intra-cavity loss

% Harmonic field parameters
Rbin = 0.7;
Rbout = 1;  %Assume no outcoupling in other mirrors
Lossb = 0.046; %Estimate of intracavity loss


%Derived quantities
tau = L./c; % Cavity round trip time

ka_in = (1-sqrt(Rain))./(tau); %Front coupler decay rate
ka_out = (1-sqrt(Raout))./(tau); %Front coupler decay rate
ka_l = (1-sqrt(1-Lossa))./(tau); %Front coupler decay rate
ka_total = ka_in + ka_out + ka_l; %The total decay rate of the whole cavity


kb_in = (1-sqrt(Rbin))./(tau); %Front coupler decay rate
kb_out = (1-sqrt(Rbout))./(tau); %Front coupler decay rate
kb_l = (1-sqrt(1-Lossb))./(tau); %Front coupler decay rate
kb_total = kb_in + kb_out + kb_l; %The total decay rate of the whole cavity

% Output some useful cavity values
% FSR = 1/tau;
% FundFinesse = (pi.*(Rain.*Raout.*(1-Lossa))^0.25)./(1-(Rain.*Raout.*(1-Lossa))^0.5);
% FundLineWidth = FSR./Finesse;

%Crystal Parameters
L_c = 10e-3; % [m]Crystal length (This may need to be tweeked to get the effective length of a gaussian focused beam through a crystal
epsilon_0 = 1090.5; % [s^-1] Non-linear coupling strength
nKTP_a = 1.830; % RI of fundamental
nKTP_b = 1.889; %RI of harmonic 
dndT_a = 1.4774*10^-5; % [1/K]Temp gradient of fundamental
dndT_b = 2.4188*10^-5;% [1/K]Temp gradient of harmonic
alpha_KTP = 6.7e-6; % [m/(m.K)] First order expation rate of KTP

%%%%%%%%%%%%%%%%%%%%%%%

% Speificies of input power etc
Deltaa_ss = 0; % On resonance no fluctations in cavity length
Deltab_ss = 0;
Deltaa_delta = 0;
Deltab_delta = 0;


epsilon_ss = epsilon_0;
epsilon_delta = 0;

a_ss = 0;
Bin = sqrt(105e-3/(h*c/lambdaHarm)); %[W] Incident pump amplitude put in units of sqrt(photon/s)
b_ss = sqrt(2.*kb_in)./(kb_total+1i.*Deltab_ss).*Bin;


%% Scan frequency as a parameter
Omega = logspace(0,10,1000);




%% Basic working model
LambdaMat = [1 1 0 0;1i -1i 0 0;0 0 1 1; 0 0 1i -1i]; %Defining transfer matrix that converts the creation anihilation form of 4x1 Fourier domain matrixes to 4x1 quadtrature terms



% Intracavity field round trip transfer matrix
Mc = [-ka_total-1i.*Deltaa_ss,epsilon_ss.*b_ss,epsilon_ss.*conj(a_ss),0;...
      conj(epsilon_ss).*conj(b_ss),-ka_total+1i.*Deltaa_ss,0,conj(epsilon_ss).*a_ss;...
      -conj(epsilon_ss).*a_ss,0,-kb_total-1i.*Deltab_ss,0;...
      0,-epsilon_ss.*conj(a_ss),0,-kb_total+1i.*Deltab_ss]; % Steady state transfer components of dual cavity equation matrix
  
% The individual mirror coupling rate matrixes
Min = diag([sqrt(2.*ka_in) sqrt(2.*ka_in) sqrt(2.*kb_in) sqrt(2.*kb_in)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Mout = diag([sqrt(2.*ka_out) sqrt(2.*ka_out) sqrt(2.*kb_out) sqrt(2.*kb_out)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Ml = diag([sqrt(2.*ka_l) sqrt(2.*ka_l) sqrt(2.*kb_l) sqrt(2.*kb_l)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths


% Input fields in quadrature form

Xin = [1;1;1;1]; % Template of what to put in [Ain_delta,AinDagger_delta,Bin_delta,BinDagger_delta]; % Input coupler fields
Xout = [1;1;1;1]; % Template of what to put in [Aout_delta,AoutDagger_delta,Bout_delta,BoutDagger_delta]; % Output coupler fields
Xl = [1;1;1;1]; % Template of what to put in [Al_delta,AlDagger_delta,Bl_delta,BlDagger_delta]; % Intracavity loss 'port' this we treat as a mirror at one point

XDelta = [-1i.*a_ss.*Deltaa_delta;1i.*a_ss.*Deltaa_delta;-1i.*b_ss.*Deltab_delta;-1i.*b_ss.*Deltab_delta]; %Flucating part of the cavity detuning, due to cavity length noise or other temperature induced detuning etc

Xepsilon = [conj(a_ss).*b_ss.*epsilon_delta;a_ss.*conj(b_ss).*conj(epsilon_delta);-0.5.*a_ss.^2.*conj(epsilon_delta);-0.5.*conj(a_ss).^2.*epsilon_delta]; %Fluctuation component of the coupling coefficient


%Compute circulating fields (this is all lumped together
%LATER THIS WILL BE SEPERATED OUT

for i = 1:length(Omega)
Xc(:,i) = inv(1i.*Omega(i).*eye(4)-Mc)*(Min*inv(LambdaMat)*Xin+Mout*inv(LambdaMat)*Xout+Ml*inv(LambdaMat)*Xl+eye(4)*XDelta+eye(4)*Xepsilon); %Fourier domain operator solved value of the intracavity circulating fields in terms of the extracavity fields
end


%Compute outcoupled fields (except for loss port as this is not accessable physically.
Xtrans = Mout*Xc-inv(LambdaMat)*Xout*ones(1,length(Xc));%Applying cavity boundtry conditions to get the transmitted and reflected fields.
Xrefl = Min*Xc-inv(LambdaMat)*Xin*ones(1,length(Xc));



QuadtratureFields_trans = LambdaMat*Xtrans;
QuadtratureFields_refl = LambdaMat*Xrefl;

V_trans = abs(QuadtratureFields_trans).^2;
V_refl = abs(QuadtratureFields_refl).^2;



figure(1)
NP = semilogx(Omega/2/pi,10*log10(V_refl(1,:)),Omega/2/pi,10*log10(V_refl(2,:)),'--');
NP(1).LineWidth = 2;
axis([min(Omega/2/pi),max(Omega/2/pi),-20,30])
legend('V1','V2')
xlabel('Frequency from Resonance [MHz]')
ylabel('V_{relf} [dBm rel SN]')
set(gca,'FontSize',18)


