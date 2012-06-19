% COSY pulse sequence.
%
% The following parameters are currently accepted:
%
% parameters.sweep              sweep width in Hz
% parameters.offset             receiver offset in Hz
% parameters.npoints_f1         number of points to be computed in f1
% parameters.npoints_f2         number of points to be computed in f2
% parameters.spins              nuclei on which the sequence runs ('1H', '13C' etc)
% parameters.angle              allows COSY45 and COSY60, radians
%
% If a Liouvillian is supplied as the third parameter, it will be used
% without questions, otherwise the Liouvillian is generated from the spin
% spin system and options given by the user during calls to create() and
% basis() functions.
%
% gareth.charnock@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function fid=cosy(spin_system,parameters,L)

% Show the banner
banner(spin_system,'sequence_banner');

%% Build the simulation infrastructure

% Set the defaults
if ~isfield(parameters,'angle')
    parameters.angle=pi/2;
end

% Print a report
report(spin_system,'cosy: computing COSY...');
sequence_report(spin_system,parameters);

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Assemble the Liouvillian
if nargin==3
    report(spin_system,'cosy: using the Liouvillian as supplied.');
else
    report(spin_system,'cosy: building the Liouvillian...');
    L=h_superop(spin_system)+1i*r_superop(spin_system)+1i*k_superop(spin_system);
end

% Apply the offset
if isfield(parameters,'offset')
    report(spin_system,'cosy: applying the offset...');
    L=L-offset(spin_system,parameters.spins,parameters.offset);
end

% Compute the coherent evolution timestep
timestep=1/parameters.sweep;

% Detect along the plus coherence
coil=state(spin_system,'L+',parameters.spins);

% Start from the thermal equilibrium
rho=equilibrium(spin_system);

% Get the pulse operators
Lp=operator(spin_system,'L+',parameters.spins);

%% Run the simulation

% Apply the first pulse
rho=step(spin_system,(Lp+Lp')/2,rho,pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints_f1-1,'trajectory');

% Flip the frequencies
rho_stack=conj(rho_stack);

% Apply the second pulse
rho_stack=step(spin_system,(Lp+Lp')/2,rho_stack,parameters.angle)-step(spin_system,(Lp-Lp')/2i,rho_stack,parameters.angle);

% Run the F2 evolution
fid=evolution(spin_system,L,coil,rho_stack,timestep,parameters.npoints_f2-1,'observable');

end

% "Bring me into the company of those who seek truth, and deliver me from
% those who have found it."
%
% Arthur C. Clarke


