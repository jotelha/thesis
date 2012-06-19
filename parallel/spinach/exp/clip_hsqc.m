% CLIP-HSQC pulse sequence.
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
% The sequence supports RDC mode -- it will use any order
% matrix supplied by the user to compute
% the HSQC including residual anisotropies of all interactions. Note that
% RDC mode is not compatible with Redfield relaxation theory.
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function [fid_P,fid_N]=clip_hsqc(spin_system,parameters,L)

%% Build the simulation infrastructure

% Print a report
report(spin_system,'clip_hsqc: computing CLIP-HSQC...');
sequence_report(spin_system,parameters);

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Assemble the Liouvillian
if exist('L','var')
    report(spin_system,'clip_hsqc: using the Liouvillian as supplied.');
else
    report(spin_system,'clip_hsqc: building the Liouvillian...');
    
    % Allow for residual anisotropy
    spin_system=residual(spin_system);
    
    [Iso,Q]=h_superop(spin_system);
    H=Iso+orientation(Q,[0 0 0]);
    
    R=r_superop(spin_system);
    K=k_superop(spin_system);
    
    L=H+1i*R+1i*K;
    Lb=H-1i*R-1i*K;
end

% Apply the offsets
if isfield(parameters,'offset_f1')&&(parameters.offset_f1~=0);
    report(spin_system,'clip_hsqc: applying the F1 offset...');
    L=L-offset(spin_system,parameters.spins_f1,parameters.offset_f1);
    Lb=Lb-offset(spin_system,parameters.spins_f1,parameters.offset_f1);
end
if isfield(parameters,'offset_f2')&&(parameters.offset_f2~=0)
    report(spin_system,'clip_hsqc: applying the F2 offset...');
    L=L-offset(spin_system,parameters.spins_f2,parameters.offset_f2);
    Lb=Lb-offset(spin_system,parameters.spins_f1,parameters.offset_f1);
end

% Compute the coherent evolution timesteps
timestep_f1=1/parameters.sweep_f1;
timestep_f2=1/parameters.sweep_f2;

% Compute the J(H,X) coupling evolution time
delta=1/(2*parameters.J);

% Compute appropriate timesteps and number of points
[delta_timestep,delta_npoints]=stepsize(L,delta/2);

% Get the pulse operators
report(spin_system,'clip_hsqc: generating the pulse operators...')
Lp_H=operator(spin_system,'L+',parameters.spins_f2);
Lm_H=operator(spin_system,'L-',parameters.spins_f2);
Lx_H=(Lp_H+Lm_H)/2; Ly_H=(Lp_H-Lm_H)/2i;
Lp_X=operator(spin_system,'L+',parameters.spins_f1);
Lm_X=operator(spin_system,'L-',parameters.spins_f1);
Lx_X=(Lp_X+Lm_X)/2;

% Get the detection state, the screening state, and the starting state
coil=state(spin_system,'L-',parameters.spins_f2);
rho=equilibrium(spin_system);

%% Run the simulation
report(spin_system,'clip_hsqc: beginning the simulation...')

% Apply the pi/2 pulse on H
rho=step(spin_system,Lx_H,rho,pi/2);

% Run the delta evolution
rho=evolution(spin_system,L,[],rho,delta_timestep,delta_npoints,'final');

% Apply the pi pulse on H and X
rho=step(spin_system,Lx_H+Lx_X,rho,pi);

% Run the delta evolution
rho=evolution(spin_system,L,[],rho,delta_timestep,delta_npoints,'final');

% Apply the pi/2 pulse on H
rho=step(spin_system,Ly_H,rho,pi/2);

% Apply the pi/2 pulse on X with phase cycling
rho_px=step(spin_system,+Lx_X,rho,pi/2);
rho_mx=step(spin_system,-Lx_X,rho,pi/2);

% Run the F1 evolution (with DSS)
destination=double(coherence(spin_system,[],[+1,-1],parameters.spins_f1)&coherence(spin_system,[],0,parameters.spins_f2));
rho_stack_px=evolution(spin_system,L,[],rho_px,timestep_f1/2,parameters.npoints_f1-1,'trajectory',destination);
rho_stack_mx=evolution(spin_system,L,[],rho_mx,timestep_f1/2,parameters.npoints_f1-1,'trajectory',destination);
rho_stack_px=step(spin_system,Lx_H,rho_stack_px,pi);
rho_stack_mx=step(spin_system,Lx_H,rho_stack_mx,pi);
rho_stack_px=evolution(spin_system,L,[],rho_stack_px,timestep_f1/2,[],'refocus',destination);
rho_stack_mx=evolution(spin_system,L,[],rho_stack_mx,timestep_f1/2,[],'refocus',destination);

% First gradient (echo-anti echo detection)
rho_stack_px=coherence(spin_system,rho_stack_px,0,parameters.spins_f2);
rho_stack_mx=coherence(spin_system,rho_stack_mx,0,parameters.spins_f2);
rho_stack_px_P=coherence(spin_system,rho_stack_px,+1,parameters.spins_f1);
rho_stack_mx_P=coherence(spin_system,rho_stack_mx,+1,parameters.spins_f1);
rho_stack_px_N=coherence(spin_system,rho_stack_px,-1,parameters.spins_f1);
rho_stack_mx_N=coherence(spin_system,rho_stack_mx,-1,parameters.spins_f1);

% Detection state run backwards through simulation
coil_stack=evolution(spin_system,Lb,[],coil,-timestep_f2,parameters.npoints_f2-1,'trajectory');

% Apply the pi/2 pulse about Lx_X with phase cycling
coil_stack_px=step(spin_system,-Lx_X,coil_stack,pi/2);
coil_stack_mx=step(spin_system,+Lx_X,coil_stack,pi/2);

% Run the delta evolution (with DSS)
destination=double(coherence(spin_system,[],-1,parameters.spins_f2)&coherence(spin_system,[],0,parameters.spins_f1));
coil_stack_px=evolution(spin_system,Lb,[],coil_stack_px,-delta_timestep,delta_npoints,'final',destination);
coil_stack_mx=evolution(spin_system,Lb,[],coil_stack_mx,-delta_timestep,delta_npoints,'final',destination);

% Second gradient
coil_stack_px=coherence(spin_system,coil_stack_px,0,parameters.spins_f1);
coil_stack_mx=coherence(spin_system,coil_stack_mx,0,parameters.spins_f1);
coil_stack_px=coherence(spin_system,coil_stack_px,-1,parameters.spins_f2);
coil_stack_mx=coherence(spin_system,coil_stack_mx,-1,parameters.spins_f2);

% Apply the pi pulse about Ly_H and Lx_X
coil_stack_px=step(spin_system,Ly_H+Lx_X,coil_stack_px,pi);
coil_stack_mx=step(spin_system,Ly_H+Lx_X,coil_stack_mx,pi);

% Run the delta evolution (with DSS)
coil_stack_px=evolution(spin_system,Lb,[],coil_stack_px,-delta_timestep,delta_npoints,'final');
coil_stack_mx=evolution(spin_system,Lb,[],coil_stack_mx,-delta_timestep,delta_npoints,'final');

% Apply the pi/2 pulse on X with phase cycling
coil_stack_px=+step(spin_system,+Lx_X,coil_stack_px,pi/2)-step(spin_system,-Lx_X,coil_stack_px,pi/2);
coil_stack_mx=+step(spin_system,-Lx_X,coil_stack_mx,pi/2)-step(spin_system,+Lx_X,coil_stack_mx,pi/2);

% Apply the pi/2 pulse on H
coil_stack_px=step(spin_system,-Lx_H,coil_stack_px,pi/2);
coil_stack_mx=step(spin_system,-Lx_H,coil_stack_mx,pi/2);

% Detection
fid_P=coil_stack_px'*rho_stack_px_P+coil_stack_mx'*rho_stack_mx_P;
fid_N=coil_stack_px'*rho_stack_px_N+coil_stack_mx'*rho_stack_mx_N;

end
