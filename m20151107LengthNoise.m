function [delta_a,delta_b] = cavityLengthNoise(omega,length)

%A function for length noise based on Kirk's and Sheila's thesis. All
%contributions of length noise come out in the fluctuating component of
%delta_a and delta_b - still unsure if static component should be included

%l = .345;%cavity length
%kb_tot = kb_in+kb_out+kb_l; %total decay rate in b

%omega = logspace(-3,6,1000);
lengthNoise = zeros(1,length(omega));

lengthNoise(1:find(omega/2/pi==100)) = 10^-14.*omega(1:find(omega/2/pi==100)).^(-1/2); %at 100Hz lengthNoise = 10^-15 m/sqrtHz
lengthNoise(find(omega/2/pi==100)+1:end)=10^-15;

lengthNoiseRMS = rms(omega',lengthNoise'); %use general spectral density/time series rms, from high to low frequency

delta_b(1) = 0; %cavity locked on fundamental, but is this different? Sheila suggests it is
delta_b(2) = omega.*lengthNoiseRMS/length; %fluctuating component of detuning
delta_a(1) = 0;
delta_a(2) = omega.*lengthNoiseRMS/length/2


%b = sqrt(2*kb_out)*Bin/(kb_tot)/(1-2*1i*delta_b/kb_tot); %Dwyer equation 5.5
% 
% THETA_B?
% 
% Mc_a = [-ka_total-1i.*Deltaa_ss, epsilon_ss.*b_ss*exp(1i*THETA_b)/(1-2*1i*delta_b/kb_tot);...
%       conj(epsilon_ss).*conj(b_ss)*exp(-1i*THETA_B)/(1+2*1i*delta_b/kb_tot),-ka_total+1i.*Deltaa_ss];...
%   %Unsure about use of conj(epsilon_ss).*conj(b_ss)

  

%calculate b
%calculate delta theta?
% 
% figure(1)
% loglog(omega,lengthNoise,omega,lengthNoiseRMS)
% 
% figure(2)
% loglog(omega,delta_b)