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
Lossa = 1-0.989; %Intra-cavity loss

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
Delta_a0 = [0 0]; % On resonance with no fluctuations in harmonic and fundamental
Delta_b0 = [0 0];

epsilon = [epsilon_0 0];

Ain = 0;
Bin = sqrt(105e-3/(h*c/lambdaHarm)); %[W] Incident pump amplitude put in units of sqrt(photon/s)



% Input fields in quadrature form
Vin = [1;1;1;1]; % Template of what to put in [Ain_delta,AinDagger_delta,Bin_delta,BinDagger_delta]; % Input coupler fields
Vout = [1;1;1;1]; % Template of what to put in [Aout_delta,AoutDagger_delta,Bout_delta,BoutDagger_delta]; % Output coupler fields


%% Example scan frequency as a parameter
Omega = logspace(0,10,1000);
% 
% [Delta_a,Delta_b] = cavityLengthNoise(Omega,L,10^-15); %trying out my shitty length noise stuff
% 
% for i = 1:length(Omega)  %with length noise
% [Vrefl1_LN(i),Vrefl2_LN(i)] = VreflTransfer(Omega(i),Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a,Delta_b,Vin,Vout);
% end
% 
% for i = 1:length(Omega)
% [Vrefl1(i),Vrefl2(i)] = VreflTransfer(Omega(i),Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a0,Delta_b0,Vin,Vout);
% end
% 
% figure(1)
% NP = semilogx(Omega/2/pi,10*log10(Vrefl1),Omega/2/pi,10*log10(Vrefl2),'--');
% hold on
% LN = semilogx(Omega/2/pi,10*log10(Vrefl1_LN),Omega/2/pi,10*log10(Vrefl2_LN),'--');
% %NP(1).LineWidth = 2;
% axis([min(Omega/2/pi),max(Omega/2/pi),-30,30])
% legend('V1','V2','V1 length noise','V2 length noise')
% xlabel('Frequency from Resonance [MHz]')
% ylabel('V_{relf} [dBm rel SN]')
% set(gca,'FontSize',18)

%% Example length noise parameter
lengthnoise = logspace(-16,1);
%Omega = 0;


for i = 1:length(lengthnoise)  %with length noise
    [Delta_a_PSD,Delta_b_PSD] = cavityLengthNoise(Omega,L,lengthnoise(i)); %trying out my shitty length noise stuff
    Delta_a(:,i) = Delta_a_PSD(:,100);
    Delta_b(:,i) = Delta_b_PSD(:,100);
    [Vrefl1_LN(i),Vrefl2_LN(i)] = VreflTransfer(Omega(100),Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a(:,i).',Delta_b(:,i).',Vin,Vout);

end

for ii = 1:length(lengthnoise)
    [Delta_a_PSD,Delta_b_PSD] = cavityLengthNoise(Omega,L,lengthnoise(ii));
    M_a_ss(ii,:) = Delta_a_PSD(1,:);
    M_a_delta(ii,:) = Delta_a_PSD(2,:);
    M_b_ss(ii,:) = Delta_b_PSD(1,:);
    M_b_delta(ii,:) = Delta_b_PSD(2,:);
end
% 
% for i = 1:length(Omega)
% [Vrefl1(i),Vrefl2(i)] = VreflTransfer(Omega(i),Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a0,Delta_b0,Vin,Vout);
% end

figure(2)
% NP = semilogx(lengthnoise,10*log10(Vrefl1),lengthnoise,10*log10(Vrefl2),'--');

LN = semilogx(lengthnoise,10*log10(Vrefl1_LN),lengthnoise,10*log10(Vrefl2_LN),'--');
%NP(1).LineWidth = 2;
axis([min(lengthnoise),max(lengthnoise),-30,30])
legend('V1','V2','V1 length noise','V2 length noise')
xlabel('Length Noise m/sqrt(Hz)')
ylabel('V_{relf} [dBm rel SN]')
set(gca,'FontSize',18)

figure(3)
semilogx(Omega/2/pi,M_a_delta, Omega/2/pi,M_b_delta)
xlabel('freq')
ylabel('a_delta')

% %% Another example, effect of changing loss
% Omega = 0;
% 
% Lossa1 = linspace(0,1,1000)
% ka_l1 = (1-sqrt(1-Lossa1))./(tau);
% 
% for ii = 1:length(Lossa1)
%     [V1refl1(ii),V1refl2(ii)] = VreflTransfer(Omega,Ain,Bin,epsilon,ka_in,ka_out,ka_l1(ii),kb_in,kb_out,kb_l,Delta_a,Delta_b,Vin,Vout);
% end
% 
% figure(2) 
% Loss = semilogx(Lossa1,10*log10(V1refl1), Lossa1,10*log10(V1refl2))
% legend('V1','V2')
% xlabel('Intracavity loss at fundamental')
% ylabel('V_{relf} [dBm rel SN]')
% set(gca,'FontSize',18)