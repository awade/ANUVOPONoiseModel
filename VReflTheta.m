function [Vrefl1,Vrefl2] = VReflTheta(THETA_in,THETA_out,THETA_loss,THETA_Delta,THETA_epsilon,Vin,Vout)
%
% Computes the 
% VReflTheta calculates Vrefl 1 and 2 from Theta matrices and Vin/Vout
% [Vtrans1,Vtrans2] = VReflTheta(THETA_in,THETA_out,THETA_loss,THETA_Delta,Theta_epsilon,Vin,Vout)
%
% Function only accomidates for injection of a seed from the same port as
% the pump. Computing transmitted field may be achieved with swapping ports
% in the input to the function.
%
% Output:
%   [Vrefl1,Vrefl2] = The quadrature field pairs relative to the pump
%   field
%
% Input: 
%   Omega = Fourier sideband frequency 
%   Vin = input port fluctuating fields
%   Vout = output port fluctuating fields
%   Vl = loss port fluctuating fields 
%
%
%
%   Ain = classical SI units power of seed at input [W]
%   Bin = classical SI units power of harmonic (pump) at input [W]
%   Delta_a = [Delta_a,Delta_a_delta] Steady state and fluctating fundmental detuning 
%   Delta_b = [Delta_a,Delta_a_delta] Steady state and fluctating harmonic detuning 
%   
%
% Notes: This was pulled out of VOPONoiseBudget_v0_1 to make its own
% function
%
% Author: Andrew Wade
% Date: 6 Nov 2015
% Mods: 18/11/2015



Vl = [1;1;1;1]; % Vacuum fields coupled throught the loss 'port'


Voutput = abs(THETA_in).^2*Vin+abs(THETA_out).^2*Vout+abs(THETA_loss).^2*Vl+abs(THETA_Delta).^2+abs(THETA_epsilon).^2; %Compute the total noise coupled through each port back through the input coupler


Vrefl1 = Voutput(1,:); Vrefl2 = Voutput(2,:); % Take the two quadrature variances and drop all the pump field flucating terms as they are uninteresting to us.




