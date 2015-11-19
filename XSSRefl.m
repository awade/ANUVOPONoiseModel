function [delta_X_ss] = XSSRefl(Omega,Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a,Delta_b,delta_X_in,delta_X_out,delta_X_l,varargin)
%
% Computes the 
% XSRefl Computes the steady state part of the quadrature delta_X_refl,
% this has been separated from the THETA_DELTA and THETA_epsilon funtions,
% where noise sources we're interested in will couple in.
%
% 
%
% Function only accomidates for injection of a seed from the same port as
% the pump. Computing transmitted field may be achieved with swapping ports
% in the input to the function.
%
% Output:
%   [delta_X_ss] = The quadrature field pairs relative to the pump
%   field
%
% Input: 
%   Omega = Fourier sideband frequency 
%   delta_X_in = input port quadtrature fluctatuons = LAMBDA*Chi_in
%   delta_X_out = output port quadtrature fluctatuons = LAMBDA*Chi_out
%   delta_X_l = loss port quadtrature fluctatuons = LAMBDA*Chi_l
%
%
%
%   Ain = classical SI units power of seed at input [W]
%   Bin = classical SI units power of harmonic (pump) at input [W]
%   Delta_a = [Delta_a,Delta_a_delta] Steady state and fluctating fundmental detuning 
%   Delta_b = [Delta_a,Delta_a_delta] Steady state and fluctating harmonic detuning 
%   
% Taken from VreflTransfer, but removing V's and other thetas and such
%
% Author: Andrew Wade
% Date: 6 Nov 2015
% Mods: 18 Nov 2015

% Check number of inputs and throw error if wrong
if nargin>14
    error('VreflTransferErr1:ErrorInputValNumberTooMany','User entered too many input variables, function requires at most 14');
end
if nargin<10
    error('VreflTransferErr2:ErrorInputValNumberTooFew','User entered too few input variables, function requires at least 10');
end
if nargin==13
    error('VreflTransferErr3:ErrorInputValNumberElevenNotRight','User entered Variance of fields at input coupler but not the output coupler, information is incompleate leave out if you want to assume they are all vacuum');
end

% Handle cases of different number of inputs
switch nargin
    case 10
        Vin = [1;1;1;1]; % Assume input port fields are vacuum
        Vout = [1;1;1;1]; % Assume out port fields are vacuum
        Delta_a = [0 0]; % Assume fundamental field is on resonance with no fluctuations
        Delta_b = [0 0]; % Assume harmonic field is on resonance with no fluctuations
    case 11
        Vin = [1;1;1;1]; % Assume input port fields are vacuum
        Vout = [1;1;1;1]; % Assume out port fields are vacuum
        Delta_b = [0 0]; % Assume harmonic field is on resonance with no fluctuations
    case 12
        Vin = [1;1;1;1]; % Assume input port fields are vacuum
        Vout = [1;1;1;1]; % Assume out port fields are vacuum
end


Vl = [1;1;1;1]; % Vacuum fields coupled throught the loss 'port'

ka_total = ka_in + ka_out + ka_l; %The total decay rate of the whole cavity
kb_total = kb_in + kb_out + kb_l; %The total decay rate of the whole cavity

% Convert variables to standard form
epsilon_ss = epsilon(1);
epsilon_delta = epsilon(2);

Deltaa_ss = Delta_a(1); % On resonance no fluctations in cavity length
Deltaa_delta = Delta_a(2);
Deltab_ss = Delta_b(1);
Deltab_delta = Delta_b(2);

b_ss = sqrt(2.*kb_in)./(kb_total+1i.*Deltab_ss).*Bin;
a_ss = sqrt(2.*ka_in).*(ka_total-1i.*Deltaa_ss+epsilon_ss.*b_ss)./(ka_total.^2+Deltaa_ss.^2-abs(epsilon_ss.*b_ss).^2).*Ain;

LambdaMat = [1 1 0 0;1i -1i 0 0;0 0 1 1; 0 0 1i -1i]; %Defining transfer matrix that converts the creation anihilation form of 4x1 Fourier domain matrixes to 4x1 quadtrature terms

% XDelta = [-1i.*a_ss.*Deltaa_delta;1i.*a_ss.*Deltaa_delta;-1i.*b_ss.*Deltab_delta;-1i.*b_ss.*Deltab_delta]; %Flucating part of the cavity detuning, due to cavity length noise or other temperature induced detuning etc
% 
% Xepsilon = [conj(a_ss).*b_ss.*epsilon_delta;a_ss.*conj(b_ss).*conj(epsilon_delta);-0.5.*a_ss.^2.*conj(epsilon_delta);-0.5.*conj(a_ss).^2.*epsilon_delta]; %Fluctuation component of the coupling coefficient


% Intracavity field round trip transfer matrix
Mc = [-ka_total-1i.*Deltaa_ss,epsilon_ss.*b_ss,epsilon_ss.*conj(a_ss),0;...
      conj(epsilon_ss).*conj(b_ss),-ka_total+1i.*Deltaa_ss,0,conj(epsilon_ss).*a_ss;...
      -conj(epsilon_ss).*a_ss,0,-kb_total-1i.*Deltab_ss,0;...
      0,-epsilon_ss.*conj(a_ss),0,-kb_total+1i.*Deltab_ss]; % Steady state transfer components of dual cavity equation matrix
  
% The individual mirror coupling rate matrixes
Min = diag([sqrt(2.*ka_in) sqrt(2.*ka_in) sqrt(2.*kb_in) sqrt(2.*kb_in)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Mout = diag([sqrt(2.*ka_out) sqrt(2.*ka_out) sqrt(2.*kb_out) sqrt(2.*kb_out)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Ml = diag([sqrt(2.*ka_l) sqrt(2.*ka_l) sqrt(2.*kb_l) sqrt(2.*kb_l)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths


%Computing the contributions to total noise sourcewise

    THETA_in = LambdaMat*(Mout*inv(1i.*Omega.*eye(4)-Mc)*Min)*inv(LambdaMat); %Compute noise coupled from different ports and disturbance mechanisms
    THETA_out = LambdaMat*(Mout*inv(1i.*Omega.*eye(4)-Mc)*Mout-eye(4))*inv(LambdaMat);
    THETA_l = LambdaMat*Mout*inv(1i.*Omega.*eye(4)-Mc)*Ml*inv(LambdaMat);
    
    delta_X_ss = THETA_in*delta_X_in + THETA_out*delta_X_out + THETA_l*delta_X_l;
    
    
    %Computing multiplication of steady state theta's and Chi's
    
    
    
    % Noting that the fields are uncorrelated they may be added seperatly
    % with no covarianace to give the ouput field variances
    %Voutput = abs(THEATA_in).^2*Vin+abs(THEATA_out).^2*Vout+abs(THEATA_l).^2*Vl+abs(THEATA_Delta).^2+abs(THEATA_epsilon).^2; %Compute the total noise coupled through each port back through the input coupler


%Vrefl1 = Voutput(1,:); Vrefl2 = Voutput(2,:); % Take the two quadrature variances and drop all the pump field flucating terms as they are uninteresting to us.



