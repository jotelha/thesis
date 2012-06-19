% The pulse-acquire sequence in zero field.
%
% The following parameters are currently accepted:
%
%       parameters.sweep - the width of the spectral window (Hz)
%       parameters.npoints  - number time steps in the simulation
%       parameters.zerofill - number of points to zerofill to
%
% If a Liouvillian L is supplied, it is used for propagation, otherwise it
% would be generated here.
%       
% ilya.kuprov@oerc.ox.ac.uk

function fid=zerofield(spin_system,parameters)

% Show the banner
banner(spin_system,'sequence_banner');

% Set the secularity assumptions
spin_system=secularity(spin_system,'lowfield');

% Assemble the Liouvillian
report(spin_system,'zerofield: building the Liouvillian...');
if strcmp(spin_system.inter.regime,'liquid')
    L=h_superop(spin_system)+1i*r_superop(spin_system)+1i*k_superop(spin_system);
else
    error('zerofield: liquid state only for now.');
end

% Magnetogyric ratio weights
weights=abs(spin_system.comp.gamma)/max(abs(spin_system.comp.gamma));

% Initialize the states and operators.
rho=zeros(spin_system.bas.nstates,1);
coil=zeros(spin_system.bas.nstates,1);
Lp=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,10*spin_system.bas.nstates);

for n=1:spin_system.comp.nspins
    
    % Initial state
    rho=rho+weights(n)*state(spin_system,'Lz',n);
    
    % Detection state
    coil=coil+weights(n)*state(spin_system,'L+',n);
    
    % Excitation operator
    Lp=Lp+weights(n)*operator(spin_system,'L+',n);
    
end

% Define Ly in terms of L+ and L-
Ly=(Lp-Lp')/2i;

% Compute the digitization parameters.
timestep=1/parameters.sweep;

% Run the simulation
rho=step(spin_system,Ly,rho,pi/2);
fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');

end






