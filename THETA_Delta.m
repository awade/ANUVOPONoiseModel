function [THETA_DELTA] = THETA_Delta(Omega,Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a,Delta_b,varargin)
%
% Computes the 
% THETA_Delta Computes the THETA_DELTA Matrix from Lambda matrix from Chi_DELTA,
% Mout, Mc and Lambda matrices. This encapsultes the time-varying
% fluctuations in the detuning in each field
%
% 
%
% Function only accomidates for injection of a seed from the same port as
% the pump. Computing transmitted field may be achieved with swapping ports
% in the input to the function.
%
% Output:
%   [THETA_DELTA] = 1x4 column vector of fluctuations in detuning. Length
%   noise goes here
%
% Input: 
%   Omega = Fourier sideband frequency 
%   Chi_Delta = 4x1 vector laid out in eq 5.17 of Kirk's thesis, 
%   Chi_Delta = [-i*a_ss*delta_DELTA_a;i*conj(a_ss)_delta_DELTA_a;...
%   ... -i*b_ss*delta_DELTA_b; i*conj(b_ss)*delta_DELTA_b]
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
%       19 Nov, added MDelta matrix and Lambda-1 to make Delta and epsilon
%       consistent with other Theta/X matrices/vectors

% Check number of inputs and throw error if wrong

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

 Chi_Delta = [-1i.*a_ss.*Deltaa_delta;1i.*a_ss.*Deltaa_delta;-1i.*b_ss.*Deltab_delta;-1i.*b_ss.*Deltab_delta]; %Flucating part of the cavity detuning, due to cavity length noise or other temperature induced detuning etc

% Intracavity field round trip transfer matrix
Mc = [-ka_total-1i.*Deltaa_ss,epsilon_ss.*b_ss,epsilon_ss.*conj(a_ss),0;...
      conj(epsilon_ss).*conj(b_ss),-ka_total+1i.*Deltaa_ss,0,conj(epsilon_ss).*a_ss;...
      -conj(epsilon_ss).*a_ss,0,-kb_total-1i.*Deltab_ss,0;...
      0,-epsilon_ss.*conj(a_ss),0,-kb_total+1i.*Deltab_ss]; % Steady state transfer components of dual cavity equation matrix
  
% The individual mirror coupling rate matrixes
Mout = diag([sqrt(2.*ka_out) sqrt(2.*ka_out) sqrt(2.*kb_out) sqrt(2.*kb_out)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Min = diag([sqrt(2.*ka_in) sqrt(2.*ka_in) sqrt(2.*kb_in) sqrt(2.*kb_in)]);   
MDelta = diag([-1i*a_ss 1i*conj(a_ss) -1i*b_ss 1i*conj(b_ss)]);  %%%% New N matrix to make Chi_delta more consistent

    THETA_DELTA = LambdaMat*Min*inv(1i.*Omega.*eye(4)-Mc)*MDelta*inv(LambdaMat); %%% Now multiplying by Lambda^-1 to turn Chi into X
    