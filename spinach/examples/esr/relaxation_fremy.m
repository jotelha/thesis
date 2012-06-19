% A Spinach transcription of EasySpin Fremy salt test file, with
% acknowledgements to Stefan Stoll.
%
% The Spinach simulation is run using explicit time propagation
% in Liouville space.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function relaxation_fremy()

% General layout
sys.magnet=0.33;
sys.regime='liquid';
sys.isotopes={'E','14N'};

% Basis set
bas.mode='ESR-2';

% Interactions
inter.zeeman.eigs=cell(2,1);
inter.zeeman.euler=cell(2,1);
inter.zeeman.eigs{1}=[2.00785 2.00590 2.00265];
inter.coupling.eigs=cell(2,2);
inter.coupling.euler=cell(2,2);
inter.coupling.eigs{1,2}=[15.4137 14.0125 80.4316]*1e6; 

% Relaxation superoperator
inter.relaxation='redfield';
inter.rlx_keep='secular';
inter.tau_c=8e-10;

% Set the sequence parameters
parameters.offset=0;
parameters.sweep=2e8;
parameters.npoints=512;
parameters.zerofill=1024;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'none-1d');

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
plot_1d(spin_system,real(spectrum),parameters);

end
