% A template file for the optimization of broadband pulses designed to be 
% resilient to the resonance mismatch and power calibration errors:
%
%    http://dx.doi.org/10.1016/j.jmr.2004.11.004
% 
% Both phases and amplitudes are varied in this template (see the "magic 
% pulse" template for the phase-only version).
%
% Optionally, the pulse may be constratined to not exceed a user-specified
% average amplitude (see the mean_amp_max parameter).
%
% A universal rotation may be obtained by setting multiple initial states
% and multiple targets (as cell arrays of state vectors).
% 
% ilya.kuprov@oerc.ox.ac.uk
% naum.gershenzon@wright.edu
%
%#ok<*PFBNS>

function broadband_pulse_template()

% Use finite-difference derivatives
sys.tols.dP_method='fd_O(h^2)';

% Set the magnetic field
sys.magnet=14.1;

% Put 100 non-interacting spins at equal intervals throughout the area
% that needs to be affected by the pulse (25 kHz either side).
n_spins=100; sys.isotopes=cell(n_spins,1);
for n=1:n_spins
    sys.isotopes{n}='13C';
end
inter.zeeman.scalar=num2cell(linspace(-166,166,n_spins));

% Select a basis set - IK-2 keeps complete basis on each spin in this case,
% but ignores multi-spin orders
bas.mode='IK-2';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalize the initial state
rho=state(spin_system,'Lz','all'); rho=rho/norm(rho);

% Set up and normalize the optimization target - this is user-selectable. The
% magnetization observable along the specified state will be maximized.
L_plus=state(spin_system,'L+','all'); L_minus=state(spin_system,'L-','all');
target=(L_plus+L_minus)/2; target=target/norm(target);

% Get the drift Liouvillian (including relaxation, if any)
L=h_superop(secularity(spin_system,'nmr'))+1i*r_superop(spin_system);

% Get the control operators
Lp=operator(spin_system,'L+','all');
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Set the power levels
n_levels=11; % Number of power levels
power_max=20e3; % Maximum power level in Hz
power_min=10e3; % Minimum power level in Hz
power_levels=2*pi*linspace(power_min,power_max,n_levels);

% Set the pulse duration (seconds) and the number of time steps
pulse_duration=1e-3;
nsteps=625;

% Set the maximum average amplitude
mean_amp_max=0.66;

% Calculate the time step
time_step=pulse_duration/nsteps;

% Take a random guess for the amplitudes (first row) and phases (second row)
guess=randn(2,nsteps);

    % Build the error functional as a sum over power levels
    function [err,grad_err]=error_function(pulse_profile)
        
        % Preallocate the arrays
        err=0; grad_err=zeros(size(pulse_profile));
        
        % Loop over the power levels
        parfor pwr=1:numel(power_levels)
            
            % Translate the phase-amplitude profile into Lx, Ly coefficients
            [x_waveform,y_waveform]=polar2cartesian(power_levels(pwr)*pulse_profile(1,:),pulse_profile(2,:));
            
            % Compute the error and gradient
            [current_err,current_grad_err]=grape(spin_system,L,{Lx,Ly},[x_waveform; y_waveform],time_step,nsteps,rho,-target);
            
            % Translate the cartesian gradient into amplitude-phase gradient
            [~,~,df_dA,df_dphi]=cartesian2polar(x_waveform,y_waveform,current_grad_err(1,:),current_grad_err(2,:));
            
            % Add to the total
            err=err+current_err; grad_err=grad_err+[power_levels(pwr)*df_dA; df_dphi];
            
        end
        
        % Normalize the objective and the gradient
        err=err/numel(power_levels); grad_err=grad_err/numel(power_levels);
        
    end

% Build the average amplitude constraint
Aneq=[ones(1,nsteps) zeros(1,nsteps)]/nsteps;

% Build the boundaries
lower_bound=[zeros(1,nsteps); -2*pi*ones(1,nsteps)];
upper_bound=[ones(1,nsteps);   2*pi*ones(1,nsteps)];

% Run the optimization
options=optimset('TolX',eps,'TolFun',eps,'Display','iter','MaxIter',500,...
                 'MaxFunEvals',Inf,'Algorithm','active-set','LargeScale','on',...
                 'GradObj','on','DerivativeCheck','off','FinDiffType','central');
answer=fmincon(@error_function,guess,Aneq,mean_amp_max,[],[],lower_bound,upper_bound,[],options);

% Plot the resulting amplitude and phase profiles
subplot(1,2,1); plot(answer(1,:));
subplot(1,2,2); plot(answer(2,:));

end

