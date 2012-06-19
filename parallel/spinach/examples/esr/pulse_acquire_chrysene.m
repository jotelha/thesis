% W-band pulsed ESR spectrum of a chrysene cation radical.
%
% ilya.kuprov@oerc.ox.ac.uk
% matthew.krzystyniak@oerc.ox.ac.uk

function pulse_acquire_chrysene()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\chrysene_cation.log'),{{'E','E'},{'H','1H'}},[],options);

sys.sym_spins={[1 7],[2 8],[3 9],[4 10],[5 11],[6 12]};
sys.sym_group={'C2','C2','C2','C2','C2','C2'};

% Set the simulation parameters
sys.magnet=3.5;
bas.mode='ESR-1';

% Set the sequence parameters
parameters.offset=-2e7;
parameters.sweep=1.5e8;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'crisp-1d');

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
plot_1d(spin_system,real(spectrum),parameters);

end

