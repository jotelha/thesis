% 27Al NMR spectrum of a static single-crystal spin system
% read in from a Simpson input file, wich acknowledgements to 
% Zdenek Tosner and Niels Chr. Nielsen.
%
% ilya.kuprov@oerc.ox.ac.uk

function simpson_test_2()

% Spin system and interactions (Simpson input)
[sys,inter]=simpson2spinach('../molecules/simpson_input_2.in');

% Extra parameters
sys.magnet=14.1;
sys.regime='crystal';
bas.mode='complete';

% Sequence parameters
parameters.offset=0;
parameters.sweep=1e7;
parameters.npoints=256;
parameters.zerofill=512;
parameters.spins='27Al';

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