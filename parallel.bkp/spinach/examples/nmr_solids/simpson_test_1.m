% 1H NMR spectrum of a static single-crystal dipolar network
% read in from a Simpson input file, wich acknowledgements to 
% Zdenek Tosner and Niels Chr. Nielsen.
%
% ilya.kuprov@oerc.ox.ac.uk

function simpson_test_1()

% Spin system and interactions (Simpson input)
[sys,inter]=simpson2spinach('../molecules/simpson_input_1.in');

% Extra parameters
sys.magnet=14.1;
sys.regime='crystal';
bas.mode='IK-1';
bas.level=6;

% Sequence parameters
parameters.offset=0;
parameters.sweep=40000;
parameters.npoints=64;
parameters.zerofill=256;
parameters.spins='1H';

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