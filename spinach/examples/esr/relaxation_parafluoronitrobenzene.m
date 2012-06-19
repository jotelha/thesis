% A Spinach transcription of EasySpin parafluoronitrobenzene test file, with
% acknowledgements to Stefan Stoll.
%
% The Spinach simulation is run using explicit time propagation
% in Liouville space.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function relaxation_parafluoronitrobenzene()

% General layout
sys.magnet=0.33898;
sys.regime='liquid';
sys.isotopes={'E','14N','19F','1H','1H','1H','1H'};
sys.sym_group={'C2','C2'};
sys.sym_spins={[4 5],[6 7]};

% Basis set
bas.mode='ESR-2';

% Zeeman interactions
inter.zeeman.eigs=cell(7,1);
inter.zeeman.euler=cell(7,1);
inter.zeeman.eigs{1}=[2.0032 2.0012 2.0097];

% Spin-spin couplings
inter.coupling.eigs=cell(7,7);
inter.coupling.euler=cell(7,7);
inter.coupling.eigs{1,2}=(40.40+[24 -12 -12])*1e6;
inter.coupling.eigs{1,3}=(22.51+[34.9 -19.8 -15])*1e6;
inter.coupling.eigs{1,4}=[9.69 9.69 9.69]*1e6;
inter.coupling.eigs{1,5}=[9.69 9.69 9.69]*1e6;
inter.coupling.eigs{1,6}=[3.16 3.16 3.16]*1e6;
inter.coupling.eigs{1,7}=[3.16 3.16 3.16]*1e6;

% Relaxation superoperator
inter.relaxation='redfield';
inter.rlx_keep='secular';
inter.tau_c=80e-12;

% Experiment parameters
parameters.offset=0;
parameters.sweep=2e8;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Perform apodization
fid=apodization(fid,'crisp-1d');

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
plot_1d(spin_system,real(spectrum),parameters);

end
