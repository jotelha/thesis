% HETCOR pulse sequence.
%
% The following parameters are currently accepted:
%
% parameters.sweep_f1           F1 sweep width in Hz
% parameters.sweep_f2           F2 sweep width in Hz
% parameters.offset_f1          F1 offset in Hz
% parameters.offset_f2          F2 offset in Hz
% parameters.npoints_f1         number of points to be computed in F1
% parameters.npoints_f2         number of points to be computed in F2
% parameters.zerofill_f1        number of points to zerofill F1 to
% parameters.zerofill_f2        number of points to zerofill F2 to
% parameters.spins_f1           heteronucleus
% parameters.spins_f2           primary nucleus
% parameters.J                  scalar coupling in Hz
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function fid=hetcor(spin_system,parameters,L)

%% Build the simulation infrastructure

% Show the banner
banner(spin_system,'sequence_banner');

% Print a report
report(spin_system,'hetcor: computing HETCOR...');
sequence_report(spin_system,parameters);

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Assemble the Liouvillian
if nargin==3
    report(spin_system,'hetcor: using the Liouvillian as supplied.');
else
    report(spin_system,'hetcor: building the Liouvillian...');
    H=h_superop(spin_system);
    R=r_superop(spin_system);
    K=k_superop(spin_system);
    L=H+1i*R+1i*K;
end

% Apply the offsets
if isfield(parameters,'offset_f1')&&(parameters.offset_f1~=0);
    report(spin_system,'hetcor: applying the F1 offset...');
    L=L-offset(spin_system,parameters.spins_f1,parameters.offset_f1);
end
if isfield(parameters,'offset_f2')&&(parameters.offset_f2~=0)
    report(spin_system,'hetcor: applying the F2 offset...');
    L=L-offset(spin_system,parameters.spins_f2,parameters.offset_f2);
end

% Compute the coherent evolution timesteps
timestep_f1=1/parameters.sweep_f1;
timestep_f2=1/parameters.sweep_f2;

% Compute the J(H,X) coupling evolution times
delta_2=1/(2*parameters.J);
delta_3=1/(3*parameters.J);

% Compute appropriate step sizes and numbers of steps
[delta_2_step, delta_2_points]=stepsize(L,delta_2);
[delta_3_step, delta_3_points]=stepsize(L,delta_3);

% Get the pulse operators
Lp_H=operator(spin_system,'L+',parameters.spins_f1);
Lm_H=operator(spin_system,'L-',parameters.spins_f1);
Lx_H=(Lp_H+Lm_H)/2; Ly_H=(Lp_H-Lm_H)/2i;
Lp_X=operator(spin_system,'L+',parameters.spins_f2);
Lm_X=operator(spin_system,'L-',parameters.spins_f2);
Lx_X=(Lp_X+Lm_X)/2;

% Get the detection state
coil=state(spin_system,'L+',parameters.spins_f2);

%% Run the simulation

% Start from the thermal equilibrium
rho=equilibrium(spin_system);

% Apply the pulse on H
rho=step(spin_system,Lx_H,rho,pi/2)+1i*step(spin_system,Ly_H,rho,pi/2);

% Run the  F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep_f1/2,parameters.npoints_f1-1,'trajectory');
rho_stack=step(spin_system,Lx_X,rho_stack,pi);
rho_stack=evolution(spin_system,L,[],rho_stack,timestep_f1/2,parameters.npoints_f1-1,'refocus');

% Run the 1/2J part of the delta evolution
rho_stack=evolution(spin_system,L,[],rho_stack,delta_2_step,delta_2_points,'final');

% Apply the pulses on H and X
rho_stack=step(spin_system,Lx_X+Lx_H,rho_stack,pi/2)+step(spin_system,Lx_X+Lx_H,rho_stack,-pi/2);

% Run the 1/3J part of the delta evolution
rho_stack=evolution(spin_system,L,[],rho_stack,delta_3_step,delta_3_points,'final');

% Decouple the H nucleus
[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.spins_f1);

% Detect
fid=evolution(spin_system,L,coil,rho_stack,timestep_f2,parameters.npoints_f2-1,'observable');

end
