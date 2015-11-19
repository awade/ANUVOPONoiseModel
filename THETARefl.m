function [THETA_in, THETA_out, THETA_loss, THETA_Delta, THETA_epsilon] = THETARefl(Omega,Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a,Delta_b,phi,varargin)
%
%
% Computes the ALL the theta matrices from 
% Mout, Mc and Lambda matrices. This encapsultes the time-varying
% fluctuations in the nonlinearity.
%
% Requires quadrature rotation function
%
% Function only accomidates for injection of a seed from the same port as
% the pump. Computing transmitted field may be achieved with swapping ports
% in the input to the function.
%
% Output:
%     [THETA_In] = 4x4 matrix of noise from input coupler
%     [THETA_Out] = 4x4 matrix of noise from output coupler
%     [THETA_loss] = 4x4 matrix of noise from loss in cavity
%     [THETA_Delta] = 1x4 column matrix of noise from input coupler
%   [THETA_epsilon] = 1x4 column vector of fluctuations in nonlinearity
%
% Input: 
%   Omega = Fourier sideband frequency 
%   
%   Ain = classical SI units power of seed at input [W]
%   Bin = classical SI units power of harmonic (pump) at input [W]
%   Delta_a = [Delta_a,Delta_a_delta] Steady state and fluctating fundmental detuning 
%   Delta_b = [Delta_a,Delta_a_delta] Steady state and fluctating harmonic detuning 
%   
% Taken from VreflTransfer, but removing V's and other thetas and such
%% Check numer of input arguments - if final element phi is missing, take phi = 0
if nargin>13
    error('THETAReflErr1:ErrorInputValNumberTooMany','User entered too many input variables, function requires at most 13');
end
if nargin<12
    error('THETAReflErr2:ErrorInputValNumberTooFew','User entered too few input variables, function requires at least 12');
end

if nargin == 12
    phi = 0;
end

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

XDelta = [-1i.*a_ss.*Deltaa_delta;1i.*a_ss.*Deltaa_delta;-1i.*b_ss.*Deltab_delta;-1i.*b_ss.*Deltab_delta]; %Flucating part of the cavity detuning, due to cavity length noise or other temperature induced detuning etc

Xepsilon = [conj(a_ss).*b_ss.*epsilon_delta;a_ss.*conj(b_ss).*conj(epsilon_delta);-0.5.*a_ss.^2.*conj(epsilon_delta);-0.5.*conj(a_ss).^2.*epsilon_delta]; %Fluctuation component of the coupling coefficient


% Intracavity field round trip transfer matrix
Mc = [-ka_total-1i.*Deltaa_ss,epsilon_ss.*b_ss,epsilon_ss.*conj(a_ss),0;...
      conj(epsilon_ss).*conj(b_ss),-ka_total+1i.*Deltaa_ss,0,conj(epsilon_ss).*a_ss;...
      -conj(epsilon_ss).*a_ss,0,-kb_total-1i.*Deltab_ss,0;...
      0,-epsilon_ss.*conj(a_ss),0,-kb_total+1i.*Deltab_ss]; % Steady state transfer components of dual cavity equation matrix
  
% The individual mirror coupling rate matrixes
Min = diag([sqrt(2.*ka_in) sqrt(2.*ka_in) sqrt(2.*kb_in) sqrt(2.*kb_in)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Mout = diag([sqrt(2.*ka_out) sqrt(2.*ka_out) sqrt(2.*kb_out) sqrt(2.*kb_out)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Ml = diag([sqrt(2.*ka_l) sqrt(2.*ka_l) sqrt(2.*kb_l) sqrt(2.*kb_l)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths


%Computing the contributions to total noise sourcewise, including rotation
%matrix

    R_phi = quadRotation(phi);

    THETA_in = R_phi*(LambdaMat*(Min*inv(1i.*Omega.*eye(4)-Mc)*Min-eye(4))*inv(LambdaMat)); %Compute noise coupled from different ports and disturbance mechanisms
    THETA_out = R_phi*(LambdaMat*Min*inv(1i.*Omega.*eye(4)-Mc)*Mout*inv(LambdaMat));
    THETA_loss = R_phi*(LambdaMat*Min*inv(1i.*Omega.*eye(4)-Mc)*Ml*inv(LambdaMat));
    THETA_Delta = R_phi*(LambdaMat*Min*inv(1i.*Omega.*eye(4)-Mc)*XDelta);
    THETA_epsilon = R_phi*(LambdaMat*Min*inv(1i.*Omega.*eye(4)-Mc)*Xepsilon);
    % Noting that the fields are uncorrelated they may be added seperatly
    % with no covarianace to give the ouput field variances
   

