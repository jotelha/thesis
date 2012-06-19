% Calculates an X axis for the Fourier-transformed spectrum based on
% the digitization parameters supplied by the user.
%
% ilya.kuprov@oerc.ox.ac.uk
% matthew.krzystyniak@oerc.ox.ac.uk

function ax=axis_1d(spin_system,parameters)

% Build the axis
ax=linspace(parameters.sweep/2,-parameters.sweep/2,parameters.zerofill);

% Apply the offset
if isfield(parameters,'offset')
    ax=ax+parameters.offset;
end

% Convert the units if necessary
switch parameters.axis_units

    case 'ppm'
        ax=1e6*(2*pi)*ax/(spin(parameters.spins)*spin_system.inter.magnet);

    case 'Gauss'
        ax=-1e4*(2*pi)*ax/spin('E')+1e4*spin_system.inter.magnet;

    case 'Hz'
        % Do nothing, we are all set already.

    otherwise
        error('axis_1d: unknown axis units.')
        
end

end

% I am now convinced that theoretical physics is actually philosophy.
%
% Max Born

