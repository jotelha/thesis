% An optimal control optimization of the excitation waveform for a powder
% spectrum of 27Al with eeQq = 4MHz quadrupole tensor.
%
% ilya.kuprov@oerc.ox.ac.uk

function oc_quadrupole()

% Set the magnetic field
sys.magnet=1.5691;

% Get the orientation grid
grid=load('..\..\exp\grids\lebedev_rank_7.dat','ASCII');
grid_size=size(grid,1); phi=pi*grid(:,1)/180; theta=pi*grid(:,2)/180;

% Set the isotopes
sys.isotopes=cell(grid_size,1);
for n=1:grid_size
    sys.isotopes{n}='27Al';
end

% Set the quadrupolar interactions (in Hz)
inter.coupling.matrix=cell(grid_size);
for n=1:grid_size
    rotation_matrix=euler2dcm(phi(n),theta(n),0);
    inter.coupling.matrix{n,n}=rotation_matrix'*[-0.1e6 0 0; 0 -0.1e6 0; 0 0 0.2e6]*rotation_matrix;
end

% Select a basis set
bas.mode='IK-2';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up the initial state
rho=equilibrium(spin_system);

% Set up the quantity to minimize
Lp=state(spin_system,'L+','all');
Lm=state(spin_system,'L-','all');
target=0.5*(Lp+Lm)/grid_size;

% Get the operators
[Iso,Q]=h_superop(secularity(spin_system,'nmr'));
L=Iso+orientation(Q,[0 0 0]);
Lp=operator(spin_system,'L+','all');
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Set the pulse parameters
max_power=250e3; %Hz
nsteps=500;
time_step=1e-7;

% Set the guess (as a fraction of maximum power)
guess=2*pi*max_power*abs(randn(2,nsteps))*1e-3;

% Compute the absolute bound on the waveform
bounds=2*pi*max_power*ones(2,nsteps);

% Set the error functional
error_function=@(waveform)grape(spin_system,L,{Lx,Ly},waveform,time_step,nsteps,rho,target);

options=optimset('TolX',0,'TolFun',0,'Display','iter','MaxIter',50,...
                 'MaxFunEvals',Inf,'Algorithm','interior-point','Hessian',{'lbfgs',50},...
                 'GradObj','on','DerivativeCheck','off','FinDiffType','central','LargeScale','off');
waveform=fmincon(error_function,guess,[],[],[],[],-bounds,bounds,[],options);

% Transform into phase-amplitude representation
[phase,amp]=cart2pol(waveform(1,:),waveform(2,:));

% Plot
figure
subplot(1,2,1); plot(180*phase/pi);
subplot(1,2,2); bar(amp);

% Set the detection state
coil=state(spin_system,'L+','all');

% Execute as a shaped pulse
rho=shaped_pulse(spin_system,L,rho,'27Al',0,180*phase/pi,amp,time_step*nsteps);

% Detection
sweep=20e6; npoints=2048; zerofill=16384;
fid=evolution(spin_system,L,coil,rho,1/sweep,npoints,'observable');

% Apodization
fid=apodization(fid,'exponential-1d',10);

% Fourier transform
spectrum=fftshift(fft(fid,zerofill));

% Plotting
figure
subplot(211), plot(-real(fid)), axis tight
subplot(212), plot(-real(spectrum)), axis tight


