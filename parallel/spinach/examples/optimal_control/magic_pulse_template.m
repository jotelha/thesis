% A template file for the "magic pulse" optimizations. The term refers to
% the family of broadband NMR pulses that are highly tolerant to resonance
% mismatch and power calibration errors:
%
%    http://dx.doi.org/10.1016/j.jmr.2005.12.010
% 
% As per the original paper, only the phases are varied. Optionally, an am-
% plitude modulation function can be set, giving a specific amplitude enve-
% lope for the pulse (a vector of ones is used by default).
%
% A universal rotation may be obtained by setting multiple initial states
% and multiple targets (as cell arrays of state vectors).
% 
% ilya.kuprov@oerc.ox.ac.uk
% naum.gershenzon@wright.edu

function magic_pulse_template()

% Use finite-difference derivatives
sys.tols.dP_method='fd_O(h^2)';

% Set the magnetic field
sys.magnet=14.1;

% Put 100 non-interacting spins at equal intervals throughout the area
% that needs to be affected by the pulse (25 kHz either side)
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

% Set the optional amplitude modulation function
amp_mod=ones(1,nsteps);

% Calculate the time step
time_step=pulse_duration/nsteps;

% Take a random guess for the phase profile
guess=randn(1,nsteps);

    % Build the error functional as a sum over power levels
    function [err,grad_err]=error_function(phase_profile)
        
        % Preallocate the arrays
        err=0; grad_err=zeros(size(phase_profile));
        
        % Loop over the power levels
        parfor pwr=1:numel(power_levels)
            
            % Translate the phase profile into Lx, Ly coefficients
            waveform=power_levels(pwr)*[amp_mod.*cos(phase_profile)
                                        amp_mod.*sin(phase_profile)];
            
            % Compute the error and gradient
            [current_err,current_grad_err]=grape(spin_system,L,{Lx,Ly},waveform,time_step,nsteps,rho,-target);
            
            % Translate the control gradient into phase gradient
            [~,~,~,current_grad_err]=cartesian2polar(waveform(1,:),waveform(2,:),current_grad_err(1,:),current_grad_err(2,:));
            
            % Add to the total
            err=err+current_err; grad_err=grad_err+current_grad_err;
            
        end
        
        % Normalize the objective and the gradient
        err=err/numel(power_levels); grad_err=grad_err/numel(power_levels);
        
    end

% Run the optimization
options=optimset('TolX',eps,'TolFun',eps,'Display','iter','MaxIter',500,...
                 'MaxFunEvals',Inf,'Algorithm','active-set','LargeScale','on',...
                 'GradObj','on','DerivativeCheck','off','FinDiffType','central');
phase_profile=fmincon(@error_function,guess,[],[],[],[],-inf(size(guess)),inf(size(guess)),[],options);

% Plot the resulting phase profile
subplot(1,2,1); plot(amp_mod.*cos(phase_profile));
subplot(1,2,2); plot(amp_mod.*sin(phase_profile));

end

