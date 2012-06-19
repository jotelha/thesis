% 1D plotting utility.
%
% ilya.kuprov@oerc.ox.ac.uk

function plot_1d(spin_system,spectrum,parameters)

% Inspect the spectrum
if nnz(imag(spectrum))>0
    report(spin_system,'plot_1d: complex spectrum received, plotting the real part...');
end

% Set the defaults
if ~isfield(parameters,'axis_units')
    if ~isfield(parameters,'nuclei')
        parameters.axis_units='Hz';
    elseif strcmp(parameters.nuclei(1),'E')
        parameters.axis_units='Gauss';
    else
        parameters.axis_units='ppm';
    end
end

% Get the axis
ax=axis_1d(spin_system,parameters);

% Compute the derivative if necessary
if isfield(parameters,'derivative')
    spectrum=fft(ifft(spectrum).*fftdiff(parameters.derivative,length(spectrum),1)');
end

% Plot the spectrum
plot(ax,real(spectrum));

% Label the axis
xlabel(parameters.axis_units);

% Invert the axis if necessary
if isfield(parameters,'invert_axis')&&parameters.invert_axis
    set(gca,'XDir','reverse');
end

end

% The only sin on earth is to do things badly.
%
% Ayn Rand

