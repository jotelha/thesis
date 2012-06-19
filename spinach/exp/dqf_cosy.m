% DQF-COSY pulse sequence.
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
%
% ilya.kuprov@oerc.ox.ac.uk
% matthew.krzystyniak@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function [fid_cos,fid_sin]=dqf_cosy(spin_system,parameters,L)

% Show the banner
banner(spin_system,'sequence_banner');

%% Build the simulation infrastructure

% Print a report
report(spin_system,'dqf_cosy: computing DQF-COSY...');
sequence_report(spin_system,parameters);

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Assemble the Liouvillian
if nargin==3
    report(spin_system,'dqf_cosy: using the Liouvillian as supplied.');
else
    report(spin_system,'dqf_cosy: building the Liouvillian...');
    L=h_superop(spin_system)+1i*r_superop(spin_system)+1i*k_superop(spin_system);
end

% Apply the offset
if isfield(parameters,'offset')
    report(spin_system,'dqf_cosy: applying the offset...');
    L=L-offset(spin_system,parameters.spins,parameters.offset);
end

% Compute the coherent evolution timestep
timestep=1/parameters.sweep;

% Get the detection state
coil=state(spin_system,'L+',parameters.spins);

% Get the pulse operators
Lp=operator(spin_system,'L+',parameters.spins);
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

%% Run the simulation

% Start from the thermal equilibrium
rho=equilibrium(spin_system);

% Apply the first pulse
rho=step(spin_system,Lx,rho,pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints_f1-1,'trajectory');

% Apply the second pulse (States hypercomplex)
rho_stack_sin=step(spin_system,Lx,rho_stack,pi/2);
rho_stack_cos=step(spin_system,Ly,rho_stack,pi/2);

% Double quantum filter
rho_stack_cos=coherence(spin_system,rho_stack_cos,[+2,-2],'all');
rho_stack_sin=coherence(spin_system,rho_stack_sin,[+2,-2],'all');

% Apply the third pulse
rho_stack_sin=step(spin_system,Lx,rho_stack_sin,pi/2);
rho_stack_cos=step(spin_system,Lx,rho_stack_cos,pi/2);

% Run the F2 evolution
fid_cos=evolution(spin_system,L,coil,rho_stack_cos,timestep,parameters.npoints_f2-1,'observable');
fid_sin=evolution(spin_system,L,coil,rho_stack_sin,timestep,parameters.npoints_f2-1,'observable');
    
end

% "For a successful technology, reality must take precedence over public
% relations, for nature cannot be fooled."
%
% Richard Feynman
