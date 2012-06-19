% Converts the [RF_amplitude, RF_phase] representation of a pulse 
% waveform and the derivatives of any function with respect to 
% those amplitudes and phases into the [RF_x, RF_y] representation
% and the derivatives of the function with respect to those X and Y
% RF values.
%
% ilya.kuprov@oerc.ox.ac.uk
% naum.gershenzon@wright.edu

function [x,y,df_dx,df_dy]=polar2cartesian(A,phi,df_dA,df_dphi)

% Transform the coordinates
x=A.*cos(phi); y=A.*sin(phi);

% Transform the derivatives
if (nargin>2)&&(nargout>2)
    df_dx=df_dA.*cos(phi)-df_dphi.*sin(phi)./A;
    df_dy=df_dA.*sin(phi)+df_dphi.*cos(phi)./A;
end

end

% A public-opinion poll is no substitute for thought.
%
% Warren Buffett

