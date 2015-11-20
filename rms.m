function rms = rms(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						%
%	find the rms value of y(x)		%
%	y(f) could be a spectral density	%
%	or y(t) could be a time series		%
%						%
%	usage: rmsY = rms(X,Y);			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = diff(squeeze(x));
dx = cat(1,[dx(1)], dx);


rms = fliplr(sqrt(cumsum(fliplr((squeeze(y).^2.*dx)'))));

