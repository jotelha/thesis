% HSQC pulse sequence.
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
% The sequence supports RDC mode -- if parameters.residual switch is set
% to 'on', it would use the order matrix supplied by the user to compute 
% the HSQC including residual anisotropies of all interactions. Note that
% RDC mode is not compatible with Redfield relaxation theory.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function [fid_P,fid_N]=hsqc(spin_system,parameters,L)

% Show the banner
banner(spin_system,'sequence_banner');

% Set the defaults
if ~isfield(parameters,'residual')
    parameters.residual='off';
end

%% Build the simulation infrastructure

% Print a report
report(spin_system,'hsqc: computing HSQC...');
sequence_report(spin_system,parameters);

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Assemble the Liouvillian
if nargin==3
    
    % Just inform the user
    report(spin_system,'hsqc: using the Liouvillian as supplied.');
    
elseif strcmp(parameters.residual,'on');
    
    % Inform the user
    report(spin_system,'hsqc: building the Liouvillian (RDC mode)...');
    
    % Apply RDC scaling
    spin_system=residual(spin_system);
    
    % Get the Liouvillian, including anisotropies
    [H_iso,RotationalBasis]=h_superop(spin_system);
    L=H_iso+orientation(RotationalBasis,[0 0 0])+...
      1i*r_superop(spin_system)+1i*k_superop(spin_system);
    
else
    
    % Inform the user
    report(spin_system,'hsqc: building the Liouvillian...');
    
    % Get the isotropic Liouvillian
    L=h_superop(spin_system)+1i*r_superop(spin_system)+1i*k_superop(spin_system);
    
end

% Apply the offsets
if isfield(parameters,'offset_f1')&&(parameters.offset_f1~=0);
    report(spin_system,'hsqc: applying the F1 offset...');
    L=L-offset(spin_system,parameters.spins_f1,parameters.offset_f1);
end
if isfield(parameters,'offset_f2')&&(parameters.offset_f2~=0)
    report(spin_system,'hsqc: applying the F2 offset...');
    L=L-offset(spin_system,parameters.spins_f2,parameters.offset_f2);
end

% Compute the coherent evolution timesteps
timestep_f1=1/parameters.sweep_f1;
timestep_f2=1/parameters.sweep_f2;

% Compute the J(H,X) coupling evolution time
delta=1/(2*parameters.J);

% Compute appropriate timesteps and number of points
[delta_timestep, delta_npoints]=stepsize(L,delta/2);

% Get the pulse operators
Lp_H=operator(spin_system,'L+',parameters.spins_f2);
Lm_H=operator(spin_system,'L-',parameters.spins_f2);
Lx_H=(Lp_H+Lm_H)/2; Ly_H=(Lp_H-Lm_H)/2i;
Lp_X=operator(spin_system,'L+',parameters.spins_f1);
Lm_X=operator(spin_system,'L-',parameters.spins_f1);
Lx_X=(Lp_X+Lm_X)/2;

% Get the starting and the detection state
coil=state(spin_system,'L+',parameters.spins_f2);
rho=equilibrium(spin_system);

%% Run the simulation

% Apply the pulse on H
rho=step(spin_system,Lx_H,rho,pi/2);

% Run the delta evolution
rho=evolution(spin_system,L,[],rho,delta_timestep,delta_npoints,'final');

% Apply the pulse on H and X
rho=step(spin_system,Lx_H+Lx_X,rho,pi);

% Run the delta evolution
rho=evolution(spin_system,L,[],rho,delta_timestep,delta_npoints,'final');

% Apply the pulse on H 
rho=step(spin_system,Ly_H,rho,pi/2);

% Apply the pulse on X with phase cycling
rho=step(spin_system,Lx_X,rho,pi/2)-step(spin_system,Lx_X,rho,-pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep_f1/2,parameters.npoints_f1-1,'trajectory');
rho_stack=step(spin_system,Lx_H,rho_stack,pi);
rho_stack=evolution(spin_system,L,[],rho_stack,timestep_f1/2,parameters.npoints_f1-1,'refocus');

% Gradient selection
rho_stack=coherence(spin_system,rho_stack,0,parameters.spins_f2);
rho_stack_P=coherence(spin_system,rho_stack,+1,parameters.spins_f1);
rho_stack_N=coherence(spin_system,rho_stack,-1,parameters.spins_f1);

% Apply the pulse on H and X
rho_stack_P=step(spin_system,Lx_H+Lx_X,rho_stack_P,pi/2);
rho_stack_N=step(spin_system,Lx_H+Lx_X,rho_stack_N,pi/2);

% Run the delta evolution
rho_stack_P=evolution(spin_system,L,[],rho_stack_P,delta_timestep,delta_npoints,'final');
rho_stack_N=evolution(spin_system,L,[],rho_stack_N,delta_timestep,delta_npoints,'final');

% Apply the pulse on H and X
rho_stack_P=step(spin_system,Lx_H+Lx_X,rho_stack_P,pi);
rho_stack_N=step(spin_system,Lx_H+Lx_X,rho_stack_N,pi);

% Run the delta evolution
rho_stack_P=evolution(spin_system,L,[],rho_stack_P,delta_timestep,delta_npoints,'final');
rho_stack_N=evolution(spin_system,L,[],rho_stack_N,delta_timestep,delta_npoints,'final');

% Decouple the X nucleus unless RDC mode is selected
if ~strcmp(parameters.residual,'on')
    [L,rho_stack_P]=decouple(spin_system,L,rho_stack_P,parameters.spins_f1);
    [L,rho_stack_N]=decouple(spin_system,L,rho_stack_N,parameters.spins_f1);
end

% Detect
fid_P=evolution(spin_system,L,coil,rho_stack_P,timestep_f2,parameters.npoints_f2-1,'observable');
fid_N=evolution(spin_system,L,coil,rho_stack_N,timestep_f2,parameters.npoints_f2-1,'observable');

end
