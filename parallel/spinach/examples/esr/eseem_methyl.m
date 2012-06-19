% X-band ESEEM spectrum of a methyl radical.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function eseem_methyl()

% Read the spin system properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g03_to_spinach(g03_parse('..\molecules\methyl.log'),{{'E','E'},{'H','1H'}},[],options);

% Set the simulation parameters
sys.magnet=0.33;
sys.regime='crystal';
sys.disable={'pt','krylov'};
bas.mode='ESR-2';

% Set the sequence parameters
parameters.npoints=512;
parameters.timestep=0.5e-8;
parameters.spins='E';
parameters.orientation=[pi/3 pi/4 pi/5];
parameters.zerofill=4096;

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
fid=eseem(spin_system,parameters);

% Plot the time domain signal
subplot(2,1,1); plot((0:(parameters.npoints-1))*parameters.timestep*1e6,real(fid)); xlabel('\mus')

% Perform crude apodization
fid=apodization(fid-mean(fid),'kaiser-1d',6);

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));
ax=linspace(-1/(parameters.timestep),1/(parameters.timestep),parameters.zerofill)*1e-6;

% Plot the spectrum
subplot(2,1,2); plot(ax,abs(spectrum)); xlabel('MHz')

end
