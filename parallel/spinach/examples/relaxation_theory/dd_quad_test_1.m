% DD-Quad cross-correlation.
%
% ilya.kuprov@oerc.ox.ac.uk

function dd_quad_test_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','14N'};

inter.coupling.scalar={0 50; 50 0};
inter.coupling.eigs{2,2}=[1e4 1e4 -2e4];
inter.coupling.euler{2,2}=[0 0 0];

inter.coordinates={[0 0 0]
                   [0 0 1.02]};

inter.relaxation='redfield';
inter.rlx_keep='secular';
inter.tau_c=1e-9;

% Basis specification
bas.mode='complete';

% Sequence parameters
parameters.sweep=500;
parameters.spins='1H';
parameters.npoints=128;
parameters.zerofill=512;
parameters.axis_units='Hz';

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