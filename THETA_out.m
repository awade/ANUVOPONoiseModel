function THETA_Out = THETA_out(Omega,Ain,Bin,epsilon,ka_in,ka_out,ka_l,kb_in,kb_out,kb_l,Delta_a,Delta_b,varargin)
%
% Computes the 
% THETA_out Computes the THETA_out Matrix from
% Mout, Mc and Lambda matrices. This encapsultes elements of the theta
% matrix to do with the output coupler
%
% Function only accomidates for injection of a seed from the same port as
% the pump. Computing transmitted field may be achieved with swapping ports
% in the input to the function.
%
% Output:
%   [THETA_Out] = 4x4 matrix

% Input: 
%   Omega = Fourier sideband frequency 

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

% Intracavity field round trip transfer matrix
Mc = [-ka_total-1i.*Deltaa_ss,epsilon_ss.*b_ss,epsilon_ss.*conj(a_ss),0;...
      conj(epsilon_ss).*conj(b_ss),-ka_total+1i.*Deltaa_ss,0,conj(epsilon_ss).*a_ss;...
      -conj(epsilon_ss).*a_ss,0,-kb_total-1i.*Deltab_ss,0;...
      0,-epsilon_ss.*conj(a_ss),0,-kb_total+1i.*Deltab_ss]; % Steady state transfer components of dual cavity equation matrix
  
% The individual mirror coupling rate matrixes
Min = diag([sqrt(2.*ka_in) sqrt(2.*ka_in) sqrt(2.*kb_in) sqrt(2.*kb_in)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
Mout = diag([sqrt(2.*ka_out) sqrt(2.*ka_out) sqrt(2.*kb_out) sqrt(2.*kb_out)]); %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
 
    THETA_Out = LambdaMat*Min*inv(1i.*Omega.*eye(4)-Mc)*Mout*inv(LambdaMat);