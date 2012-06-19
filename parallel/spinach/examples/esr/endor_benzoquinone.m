% CW ENDOR on 2-methoxy-1,4-benzoquinone radical in liquid state.
%
% H. Joela et al., Magn. Reson. Chem., 28 (1990) 261-267.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function endor_benzoquinone()

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

% Set the sequence parameters
parameters.type='cw';
parameters.sweep=50e6;
parameters.npoints=1024;
parameters.electron_spin='E';
parameters.nuclear_spins='1H';
parameters.zerofill=4096;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=endor(spin_system,parameters);

% Perform crude apodization
fid=apodization(fid-mean(fid),'kaiser-1d',20);

% Perform Fourier transform
spectrum=fftshift(fft(real(fid),parameters.zerofill));
ax=linspace(-parameters.sweep/2,parameters.sweep/2,parameters.zerofill)*1e-6;

% Plot the spectrum
plot(ax,abs(spectrum)); xlabel('MHz')

end