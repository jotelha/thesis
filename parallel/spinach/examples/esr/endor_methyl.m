% Mims ENDOR spectrum of a methyl radical in liquid state.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function endor_methyl()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\methyl.log'),{{'E','E'},{'H','1H'}},[],options);

% Set the simulation parameters
sys.magnet=0.33;
sys.sym_group={'S3'};
sys.sym_spins={[1 2 3]};
sys.regime='liquid';
bas.mode='complete';

% Set the sequence parameters
parameters.type='mims';
parameters.npoints=512;
parameters.sweep=4e8;
parameters.tau=100e-9;
parameters.zerofill=4096;
parameters.electron_spin='E';
parameters.nuclear_spins='1H';

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=endor(spin_system,parameters);

% Perform crude apodization
fid=apodization(fid-mean(fid),'kaiser-1d',6);

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));
ax=linspace(-parameters.sweep/2,parameters.sweep/2,parameters.zerofill)*1e-6;

% Plot the spectrum
plot(ax,abs(spectrum)); xlabel('MHz')

end