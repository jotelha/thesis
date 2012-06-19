% Hypercomplex NOESY pulse sequence.
%
% The following parameters are currently accepted:
%
% parameters.sweep              sweep width in Hz
% parameters.offset             offset in Hz
% parameters.npoints_f1         number of points to be computed in f1
% parameters.npoints_f2         number of points to be computed in f2
% parameters.zerofill_f1        number of points to zerofill f1 to
% parameters.zerofill_f2        number of points to zerofill f2 to
% parameters.spins              nuclei on which the sequence runs ('1H', '13C' etc)
% parameters.tmix               NOESY mixing time in seconds
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function [cos_term,sin_term]=noesy(spin_system,parameters,L)

% Show the banner
banner(spin_system,'sequence_banner');

%% Build the simulation infrastructure

% Print a report
report(spin_system,'noesy: computing NOESY...');
sequence_report(spin_system,parameters);

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Assemble the Liouvillian
if nargin==3
    report(spin_system,'noesy: using the Liouvillian as supplied.');
else
    report(spin_system,'noesy: building the Liouvillian...');
    H=h_superop(spin_system);
    R=r_superop(spin_system);
    K=k_superop(spin_system);
    L=H+1i*R+1i*K;
end

% Apply the offset
if isfield(parameters,'offset')
    report(spin_system,'noesy: applying the offset...');
    L=L-offset(spin_system,parameters.spins,parameters.offset);
end

% Compute the coherent evolution timestep
timestep=1/parameters.sweep;

% Compute appropriate tmix timestep and number of points
[tmix_step,tmix_npoints]=stepsize(R,parameters.tmix);

% Get the detection state
coil=state(spin_system,'L+',parameters.spins);

% Get the pulse operators
Lp=operator(spin_system,'L+',parameters.spins);
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

%% Run the simulation in hypercomplex mode

% Start from the thermal equilibrium
rho=equilibrium(spin_system);

% Apply the first pulse
rho=step(spin_system,Lx,rho,pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints_f1-1,'trajectory');

% Run the second pulse
rho_stack_cos=step(spin_system,Lx,rho_stack,pi/2);
rho_stack_sin=step(spin_system,Ly,rho_stack,pi/2);

% Run homospoil
rho_stack_cos=homospoil(spin_system,rho_stack_cos);
rho_stack_sin=homospoil(spin_system,rho_stack_sin);

% Run the mixing time
rho_stack_cos=evolution(spin_system,1i*R,[],rho_stack_cos,tmix_step,tmix_npoints,'final');
rho_stack_sin=evolution(spin_system,1i*R,[],rho_stack_sin,tmix_step,tmix_npoints,'final');

% Run the third pulse
rho_stack_a=step(spin_system,Ly,rho_stack_cos,pi/2);
rho_stack_b=step(spin_system,Ly,rho_stack_sin,pi/2);

% Run the F2 evolution
cos_term=evolution(spin_system,L,coil,rho_stack_a,timestep,parameters.npoints_f2-1,'observable');
sin_term=evolution(spin_system,L,coil,rho_stack_b,timestep,parameters.npoints_f2-1,'observable');

end

% According to a trade legend, Anil Kumar had to run the very first NOESY
% experiment on a Saturday -- his supervisors viewed it as a waste of
% valuable spectrometer time.

