function [delta_a,delta_b,lengthNoiseRMS,lengthNoise] = cavityLengthNoise(FFreq,lambdaFund,lambdaHarm,CavityLength,lengthNoisemRtHz)

%A function for length noise based on Kirk's and Sheila's thesis. All
%contributions of length noise come out in the fluctuating component of
%delta_a and delta_b - still unsure if static component should be included

%will put in a variable length noise soon
% 
% Omega = logspace(0,10,1000);
% L = 0.345;
% lengthNoisemRtHz = 10^-10;

lengthNoise = zeros(1,numel(FFreq));
c = 3e8; 

[c1,f100] = min(abs(FFreq-100)); %find nearest point to 100Hz for our crossover
% % %f100 = Omega(index);

lengthNoise(1:f100) = 10*lengthNoisemRtHz.*(FFreq(1:f100)).^(-1/2); %at 100Hz lengthNoise = 10^-15 m/sqrtHz
lengthNoise(f100+1:end) = lengthNoisemRtHz;

lengthNoiseRMS = rms(FFreq',lengthNoise'); %use general spectral density/time series rms, from high to low frequency

delta_b = zeros(2,numel(FFreq));
delta_a = zeros(2,numel(FFreq));

delta_b(1,:) = 0; %cavity locked on fundamental, but is this different? Sheila suggests it is
delta_b(2,:) = (2*pi*c/lambdaHarm).*lengthNoiseRMS/CavityLength; %fluctuating component of detuning
delta_a(1,:) = 0;
delta_a(2,:) = (2*pi*c/lambdaFund).*lengthNoiseRMS/CavityLength/2;
% 
% figure
% loglog(Omega/2/pi,lengthNoise)