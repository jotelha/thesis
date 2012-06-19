% Pulse-acquire sequence in a system involving a high-spin nucleus.
%
% ilya.kuprov@oerc.ox.ac.uk

function high_spin_test_1()

% General setup
sys.magnet=14.1;
bas.mode='complete';

% Spin system
sys.isotopes={'1H','235U','1H','1H'};
inter.zeeman.scalar={-0.5  0.0  2.5  1.3};
inter.coupling.scalar{1,2}=100;
inter.coupling.scalar{3,4}=50;
inter.coupling.scalar{4,4}=0;

% Pulse sequence parameters
parameters.sweep=3500; %Hz
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.spins='1H';
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'crisp-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
plot_1d(spin_system,real(spectrum),parameters);

end

