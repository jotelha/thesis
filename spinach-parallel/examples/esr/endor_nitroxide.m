% Mims ENDOR on a 15N-labelled nitroxide radical in liquid state.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function endor_nitroxide()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\nitroxide.log'),{{'E','E'},{'N','15N'}},[],options);

% Set the simulation parameters
sys.magnet=0.33;
sys.regime='liquid';
bas.mode='complete';

% Set the sequence parameters
parameters.type='mims';
parameters.npoints=512;
parameters.sweep=1e8;
parameters.tau=100e-9;
parameters.zerofill=4096;
parameters.electron_spin='E';
parameters.nuclear_spins='15N';

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
