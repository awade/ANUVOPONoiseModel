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
Lossa = 1-0;989; 0.0; %Intra-cavity loss

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





%% Basic working model

% Intracavity field round trip transfer matrix
Mc = [-ka_total-1i.*Delta_a_ss,epsilon_ss.*b_ss,epsilon_ss.*conj(a_ss),0;...
      conj(epsilon_ss).*conj(b_ss),-ka_total+1i.*Delta_a_ss,0,conj(epsilon_ss).*a_ss;...
      -conj(epsilon_ss).*a_ss,0,-kb_total-1i.*Delta_b_ss,0;...
      0,-epsilon_ss.*conj(a_ss),0,-kb_total+1i.*Delta_b_ss]; % Steady state transfer components of dual cavity equation matrix
  
% The indivisual mirror coupling rate matrixes
Min = diag([sqrt(2.*ka_in) sqrt(2.*ka_in) sqrt(2.*kb_in) sqrt(2.*kb_in)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Mout = diag([sqrt(2.*ka_out) sqrt(2.*ka_out) sqrt(2.*kb_out) sqrt(2.*kb_out)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Ml = diag([sqrt(2.*ka_l) sqrt(2.*ka_l) sqrt(2.*kb_l) sqrt(2.*kb_l)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths

  
