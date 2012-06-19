% A Spinach transcription of EasySpin biaryl test file, with
% acknowledgements to Stefan Stoll.
% 
% The Spinach simulation is run using explicit time propagation
% in Liouville space.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_biaryl()

% General layout
sys.magnet=0.33;
sys.isotopes={'E','14N','1H','1H','1H','1H','1H','14N','1H','1H','1H','1H','1H'};
sys.sym_group={'C2','C2','C2','C2','C2','C2'};
sys.sym_spins={[2 8],[3 9],[4 10],[5 11],[6 12],[7 13]};

% Basis set
bas.mode='ESR-1';

% Zeeman interactions
inter.zeeman.scalar=cell(13,1);
inter.zeeman.scalar{1}=2.00316;

% Spin-spin couplings
inter.coupling.scalar=cell(13,13);
inter.coupling.scalar{1,2} = 12.16e6;
inter.coupling.scalar{1,3} = -6.7e6;
inter.coupling.scalar{1,4} = -1.82e6;
inter.coupling.scalar{1,5} = -7.88e6;
inter.coupling.scalar{1,6} = -0.64e6;
inter.coupling.scalar{1,7} = 67.93e6;
inter.coupling.scalar{1,8} = 12.16e6;
inter.coupling.scalar{1,9} = -6.7e6;
inter.coupling.scalar{1,10} = -1.82e6;
inter.coupling.scalar{1,11} = -7.88e6;
inter.coupling.scalar{1,12} = -0.64e6;
inter.coupling.scalar{1,13} = 67.93e6;

% Experiment parameters
parameters.offset=0;
parameters.sweep=3e8;
parameters.npoints=4096;
parameters.zerofill=8192;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Spinach code
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
