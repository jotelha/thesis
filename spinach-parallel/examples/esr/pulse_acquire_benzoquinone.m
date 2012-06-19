% ESR on 2-methoxy-1,4-benzoquinone radical in liquid state.
%
% H. Joela et al., Magn. Reson. Chem., 28 (1990) 261-267.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function pulse_acquire_benzoquinone()

% Set the simulation parameters
sys.isotopes={'E','1H','1H','1H','1H','1H','1H'};
sys.sym_group={'S3'};
sys.sym_spins={[2 3 4]};
sys.magnet=0.33;
sys.regime='liquid';

% Set the interactions
inter.zeeman.scalar={2.004577 0 0 0 0 0 0};
inter.coupling.scalar=cell(7,7);
inter.coupling.scalar{2,1}=mt2hz(0.08); 
inter.coupling.scalar{3,1}=mt2hz(0.08); 
inter.coupling.scalar{4,1}=mt2hz(0.08); 
inter.coupling.scalar{5,1}=mt2hz(-0.059); 
inter.coupling.scalar{6,1}=mt2hz(-0.364); 
inter.coupling.scalar{7,1}=mt2hz(-0.204); 

% Set the basis
bas.mode='complete';

% Experiment parameters
parameters.offset=-1e7;
parameters.sweep=3e7;
parameters.npoints=1024;
parameters.zerofill=8192;
parameters.spins='E';
parameters.axis_units='Gauss';
parameters.derivative=1;

% Spinach code
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=pulse_acquire(spin_system,parameters);

% Apodization
fid=apodization(fid,'exponential-1d',20);

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
plot_1d(spin_system,real(spectrum),parameters);

end