% Optimization of selective pulses by varying the coefficients of the 
% user-specified waveform basis functions. The treatment closely follows
% the paper by Veshtort and Griffin:
%
%     http://dx.doi.org/10.1002/cphc.200400018
%
% but uses BFGS-GRAPE instead of grid search for optimization and also 
% makes the pulses resilient to the inevitable B1 field inhomogeneity.
%
% A richer library of waveform basis functions is also available.
%
% ilya.kuprov@oerc.ox.ac.uk
% naum.gershenzon@wright.edu

function selective_pulse_template()

% Use finite-difference derivatives
sys.tols.dP_method='fd_O(h^2)';

% Set the magnetic field
sys.magnet=14.1;

% Put 100 non-interacting protons at equal intervals throughout a 4 kHz area
n_spins=100; sys.isotopes=cell(n_spins,1);
for n=1:n_spins
    sys.isotopes{n}='1H';
end
inter.zeeman.scalar=num2cell(linspace(-3.33,3.33,n_spins));

% Select a basis set - IK-2 keeps complete basis on each spin in this case,
% but ignores multi-spin orders
bas.mode='IK-2';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalize the initial state
rho=state(spin_system,'Lz','all'); rho=rho/norm(rho);

% Central 1000 Hz should be excited, 150 Hz on either side unspecified, the rest left intact
target=weighted_target(spin_system,'Lx',[zeros(1,38) ones(1,24) zeros(1,38)])+...
       weighted_target(spin_system,'Lz',[ones(1,34) zeros(1,32) ones(1,34)]);
target=target/norm(target);

% Get the drift Liouvillian (including relaxation, if any)
L=h_superop(secularity(spin_system,'nmr'))+1i*r_superop(spin_system);

% Get the control operators
Lp=operator(spin_system,'L+','all');
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Set the number of basis functions
n_functions=20;

% Set the power levels
n_levels=4; % Number of power levels
power_level=0.5e3; % Nominal power level in Hz
power_range=1.2; % Power variation
power_levels=2*pi*linspace(power_level/power_range,power_level*power_range,n_levels);

% Set the pulse duration (seconds) and the number of time steps
pulse_duration=1e-2; nsteps=200;

% Calculate the time step
time_step=pulse_duration/nsteps;

% Get the basis set
harmonics=wave_basis('sine_waves',n_functions,nsteps);

% Take a random guess for the amplitudes of harmonics
guess=randn(n_functions,2)/100;

    % Build the error functional as a sum over power levels
    function [err,grad_err]=error_function(harmonic_amplitudes)
        
        % Preallocate the arrays
        err=0; grad_err=zeros(size(harmonic_amplitudes));
        
        % Loop over the power levels
        parfor pwr=1:numel(power_levels)
            
            % Translate the harmonic amplitudes into Lx, Ly coefficients
            waveform=power_levels(pwr)*(harmonics*harmonic_amplitudes)';
                                    
            % Compute the GRAPE error and error gradient
            [current_err,current_grad_err]=grape(spin_system,L,{Lx,Ly},waveform,time_step,nsteps,rho,-target);
            
            % Translate the control gradient into harmonic amplitude gradient
            current_grad_err=power_levels(pwr)*(current_grad_err*harmonics)';
            
            % Add the chunk to the total
            err=err+current_err; grad_err=grad_err+current_grad_err;
            
        end
        
        % Normalize the error and the gradient
        err=err/numel(power_levels);
        grad_err=grad_err/numel(power_levels);
        
    end

    % Build the constraint function
    function [c,ceq]=constraint_function(harmonic_amplitudes)
        
        % Get the normalized waveform modulus
        waveform=abs(harmonics*harmonic_amplitudes);
        
        % Get the violation signal
        c=max(waveform(:))-1; ceq=0;
        
    end

% Run constrained memory-conserving BFGS
options=optimset('TolX',eps,'TolFun',eps,'Display','iter','MaxIter',50,...
                 'MaxFunEvals',Inf,'Algorithm','active-set','LargeScale','on',...
                 'GradObj','on','DerivativeCheck','off','FinDiffType','central');
harmonic_amplitudes=fmincon(@error_function,guess,[],[],[],[],[],[],@constraint_function,options);

% Plot the resulting phase profile
subplot(1,2,1); plot(harmonics*harmonic_amplitudes(:,1));
subplot(1,2,2); plot(harmonics*harmonic_amplitudes(:,2));

end

