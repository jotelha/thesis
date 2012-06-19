% Converts the [RF_x, RF_y] representation of a pulse waveform and
% the derivatives of any function with respect to those RF values
% into the [RF_amplitude, RF_phase] representation and the derivatives
% of the function with respect to the amplitudes and phases.
%
% ilya.kuprov@oerc.ox.ac.uk
% naum.gershenzon@wright.edu

function [A,phi,df_dA,df_dphi]=cartesian2polar(x,y,df_dx,df_dy)

% Transform the coordinates
A=sqrt(x.^2+y.^2); phi=atan2(y,x);

% Transform the derivatives
if (nargin>2)&&(nargout>2)
    df_dA=df_dx.*cos(phi)+df_dy.*sin(phi);
    df_dphi=-df_dx.*y+df_dy.*x;
end

end

% Beware of geeks bearing formulas.
%
% Warren Buffett 

