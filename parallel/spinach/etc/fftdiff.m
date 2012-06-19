% Spectral differentiation kernel.
%
% ilya.kuprov@oerc.ox.ac.uk

function kern=fftdiff(order,npoints,dx)

kern=ifftshift((2i*pi*((1-npoints)/2:(npoints/2))/(npoints*dx)).^order);

end

% "The aim of science is to make difficult things understandable 
% in a simpler way; the aim of poetry is to state simple things 
% in an incomprehensible way. The two are incompatible."
%
% Paul Dirac