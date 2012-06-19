% Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain
% simulation in Liouville space with averaging over a Lebedev grid.
%
% Source: H.L. Flanagan, D.J Singel, J Chem Phys, 87(10) 1987 5606-5616.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function eseem_nitroxide_14n()

% Set the simulation parameters
sys.isotopes={'14N','E'};
sys.magnet=0.3249;
sys.regime='powder';
sys.tols.grid_rank=35;
sys.disable={'pt','krylov'};
bas.mode='complete';

% Set the coupling parameters
inter.coupling.eigs=cell(2,2);
inter.coupling.eigs{1,1}=[-0.4 -1.6 +2.0]*1e5;
inter.coupling.eigs{1,2}=[2.0 2.0 2.0]*1e6;
inter.coupling.euler=cell(2,2);

% Set the sequence parameters
parameters.npoints=512;
parameters.timestep=2e-7;
parameters.spins='E';
parameters.zerofill=2048;

% Generate spin_system
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Run Spinach
fid=eseem(spin_system,parameters);

% Run apodization
fid=apodization(fid-mean(fid),'exponential-1d',5);

% Run Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the results
subplot(2,1,1);
time_axis=(0:(parameters.npoints-1))*parameters.timestep*1e6;
plot(time_axis,real(fid)); xlabel('\mus')
subplot(2,1,2);
freq_axis=linspace(-1/parameters.timestep,1/parameters.timestep,parameters.zerofill)*1e-6;
plot(freq_axis,real(spectrum)); xlabel('MHz')

end
