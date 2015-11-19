function [Vrefl1,Vrefl2] = VReflTheta(varargin) %THETA_in,THETA_out,THETA_loss,THETA_Delta,THETA_epsilon,Vin,Vout)
%
% Computes the 
% VReflTheta calculates Vrefl 1 and 2 from Theta matrices and thier associated input quadrature fluctations Vin/Vout
% [Vtrans1,Vtrans2] = VReflTheta(THETA_1,X_1,...Theta_n,X_n,Phi)
%
%
% Here the function provides the varaiance from the sum of the uncorrlated
% X_i noises transfered through the cavity system via each associated
% THETA_i function.  A quadrature of measurment is selected with the basis
% angle Phi (defult zero) that that allows selection of any abitrary pair
% of orthgonal quadratures to measuere in.
%
% For reference and further reading on the topic check out Kirk McKenzie's
% thesis (2008, PhD dissertation, The Australian National University,
% Canberra, Australia).  Also similar but slightly different ac see Shelia
% Dwyer's thesis (2013, PhD dissertation, Massachusetts Institute of
% Technology, Cambridge Massachusetts, USA)
%
%
% Output:
%   [Vrefl1,Vrefl2] = The quadrature field pairs relative to the pump
%   field
%
% Input: 
%   THETA_i = 4x4 transfer matrix of quantum noise from input port to the
%   outcoupled port for fundamental and pump. Minimum of one of these must
%   be provided, else an error will be thrown. Any number may be entered
%   into the function provided X_i for each is provided.
%
%   X_i = 4x1 Input port variance for noise source for fundamental and pump
%   field. Must be an X_i for every THETA_i inputed.
%
%   Phi = Angle of quadrature measurment basis for computing variance.
%   [Rad]
%
%
%   
%
% Notes: 
%
% Authors: Andrew Wade/Geogia Mansell
% Date: 6 Nov 2015
% Mods: 18/11/2015
% Mods: 19/11/2015 AWade generaised the script to take an abitrary number
% of \Theta_transfer matrixes and their quadrature intitial noise sources
% as well as the detection angle as well. I have also tried to include some
% basic error management to throw an error message if users does a bad.


% Error checking to ensure user is not an idiot
ErrTooFewArg.message = 'Not enough input arguments.  Minimum of one THETA and accompaning noise term X must be provided'; ErrTooFewArg.identifier =  'VReflThetaError:Error1TooFewInputs';
ArgumentsWrongShape.message = 'Arguments provided do not have the correct order and/or shape, function requires pairs of 4x4 matrixes follwed by the input noise terms 4x1, an optional quadrature phase may be then be speificed as a scaler value in Rads';ArgumentsWrongShape.identifier = 'VReflThetaError:Error2WrongInputs';
% Generic New Errror... ErrName.message = 'Message to user'; ErrName.identifier = 'FunctionName:ErrorName':

if nargin <2
    error(ErrTooFewArg) % Throw error if too few arguments,minimum of two is needed
end

%Throw error if the dimentions of provided THETA_i's and X_i's are not
%correct.  
for ii = 1:floor(nargin/2) %Step through pairs of THETA_i and X_i, skipping possible non-defult value of Phi when user specifies
  if ~isequal(size(varargin{ii*2-1}),[4 4]) || ~isequal(size(varargin{ii*2}),[4 1]) % Check dimention of THETA and or X matrix is correct, otherwise throws error
  error(ArgumentsWrongShape) % Throws error back to user with explaination
  end
end


if mod(nargin,2)==0 % Defults the value of Phi to zero in case that only even number of arguments are provided, this implies only THETA and X are provided
    Phi = 0;
else 
    Phi = varargin{nargin}; % Sets phi to the user specified in case that it is provided.
end


%Loop stepping through all the inputed THETA_i and X_i providing a rotation
%and computing the variance of each noise term.  These are later summed
%accross for total noise.
Vout_Terms = zeros(4,floor(nargin/2)); % Generate and pad with zeros for loop effiecency
for ii = 1:floor(nargin/2) %Loop stepping through each pair of THETA and X terms
    Vout_Terms(:,ii) = abs((varargin{ii*2-1})*quadRotation(Phi)).^2*varargin{ii*2}; %Compute matrix with variance entries for each noise source
end

Vrefl1 = sum(Vout_Terms(1,:));%Select first quadrature of total output noise for the fundamental and sum accross all the noises computed above.
Vrefl2 = sum(Vout_Terms(1,:));%Select second quadrature of total output noise for the fundamental and sum accross all the noises computed above.

