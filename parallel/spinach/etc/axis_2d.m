% Calculates the F1 and F2 axes for the Fourier-transformed 2D NMR spectrum 
% based on the digitization parameters supplied by the user. Syntax:
%
%   [a_f1,a_f2]=axis_2d(spin_system,parameters,center_offset,units)
%
% where parameters are the same as those supplied to pulse_acquire.m 
% and units are 'ppm', 'Hz' or 'rad/s'.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function [a_f1,a_f2]=axis_2d(spin_system,parameters)

% Set defaults
if ~isfield(parameters,'offset_f1')
    parameters.offset_f1=0;
end
if ~isfield(parameters,'offset_f2')
    parameters.offset_f2=0;
end

% Accommodate homonuclear experiments
if isfield(parameters,'offset')
    parameters.offset_f1=parameters.offset;
    parameters.offset_f2=parameters.offset;
end
if isfield(parameters,'sweep')
    parameters.sweep_f1=parameters.sweep;
    parameters.sweep_f2=parameters.sweep;
end
if isfield(parameters,'spins')
    parameters.spins_f1=parameters.spins;
    parameters.spins_f2=parameters.spins;
end

% Build the axes
a_f1=linspace(parameters.sweep_f1/2,-parameters.sweep_f1/2,parameters.zerofill_f1)+parameters.offset_f1;
a_f2=linspace(parameters.sweep_f2/2,-parameters.sweep_f2/2,parameters.zerofill_f2)+parameters.offset_f2;

% Convert the units
switch parameters.axis_units
    case 'ppm'
        % Pull the carrier frequency from the data structure
        gamma_1=spin(parameters.spins_f1);
        a_f1=1e6*(2*pi)*a_f1/(gamma_1*spin_system.inter.magnet);
        gamma_2=spin(parameters.spins_f2);
        a_f2=1e6*(2*pi)*a_f2/(gamma_2*spin_system.inter.magnet);
    case 'Hz'
        % Do nothing, we are all set already.
    case 'rad/s'
        a_f1=2*pi*a_f1;
        a_f2=2*pi*a_f2;
    otherwise
        error('axis_2d: unknown axis units');
end

end

% "Studies have suggested that hypomania can heighten certain cognitive
% processes, increase original and idiosyncratic thought, and even enhance
% linguistic skills. Manic states can also reduce the need for sleep,
% foster intense and obsessive concentration, create unmitigated
% self-confidence, and eliminate concern for social norms - just what you
% need, perhaps, to push the envelope of artistic creativity."
%
% M.F. Bear, B.W. Connors, M.A. Paradiso,
% "Neuroscience: exploring the brain", Kluwer, 2001.