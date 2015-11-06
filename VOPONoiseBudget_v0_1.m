%% NOISE BUDGET SCRIPT FOR ANU VOPO PROJECT

% Authours: Andrew Wade
% Date:
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
h = 6.626e-34 %[J.s]


Mc = [-ka-1i.*Delta_a_ss,epsilon_ss.*b_ss,epsilon_ss.*conj(a_ss),0;...
      conj(epsilon_ss).*conj(b_ss),-ka+1i.*Delta_a_ss,0,conj(epsilon_ss).*a_ss;...
      -conj(epsilon_ss).*a_ss,0,-kb-1i.*Delta_b_ss,0;...
      0,-epsilon_ss.*conj(a_ss),0,-kb+1i.*Delta_b_ss]; % Steady state transfer components of dual cavity equation matrix
  
Mj = diag([]) %Generate diagonal matrix with coupling rates of port of OPO cavity for both wavelengths
  
  
