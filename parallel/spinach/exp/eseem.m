% ESEEM pulse sequence.
%
% The following parameters are currently accepted:
%
% parameters.npoints            number of points to be computed
% parameters.timestep           simulation time step
% parameters.spins              working electron location
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function fid=eseem(spin_system,parameters)

% Show the banner
banner(spin_system,'sequence_banner');

%% Build the simulation infrastructure

% Print a report
report(spin_system,'eseem: computing ESEEM...');
sequence_report(spin_system,parameters);

% Get the Liouvilian
report(spin_system,'eseem: building the Liouvillian...');
[Iso,Q]=h_superop(secularity(spin_system,'eseem'));
R=r_superop(spin_system);

switch spin_system.inter.regime

    case 'crystal'

        H=Iso+orientation(Q,parameters.orientation);
        K=k_superop(spin_system);
        L=H+1i*R+1i*K;

        % Get the pulse operators
        Lp=operator(spin_system,'L+',parameters.spins);
        Lm=operator(spin_system,'L-',parameters.spins);
        Ly=(Lp-Lm)/2i;

        % Get the source, the detection state, and the screening state
        rho=equilibrium(spin_system);
        coil=state(spin_system,'L+',parameters.spins);
        screen=state(spin_system,'L-',parameters.spins);
        
    case 'powder'
        
        % Get the spherical averaging grid
        grid=load([spin_system.sys.root_dir spin_system.sys.slash 'exp' ...
                                            spin_system.sys.slash 'grids' ...
                                            spin_system.sys.slash 'lebedev_rank_' num2str(spin_system.tols.grid_rank) '.dat'],'ASCII');
        grid_size=size(grid,1); phi=pi*grid(:,1)/180; theta=pi*grid(:,2)/180; weight=grid(:,3);
        
        % Get the orientation array
        L_aniso=orientation(Q,[phi theta zeros(size(theta))]);
        L=blkdiag(L_aniso{:})+kron(speye(grid_size),Iso);
        L=clean_up(spin_system,L,spin_system.tols.liouv_zero);
        
        % Get the pulse operators
        Lp=kron(speye(grid_size),operator(spin_system,'L+',parameters.spins));
        Lm=kron(speye(grid_size),operator(spin_system,'L-',parameters.spins));
        Ly=(Lp-Lm)/2i;
        
        % Get the initial, screen, and detection states
        rho=kron(ones(grid_size,1),equilibrium(spin_system));
        screen=kron(weight,state(spin_system,'L-',parameters.spins));
        coil=kron(weight,state(spin_system,'L+',parameters.spins));
        
    otherwise
        
        error('eseem: sys.regime not set');
        
end

%% Run the simulation

% Apply the first pulse
report(spin_system,'eseem: applying the first pulse...')
rho=step(spin_system,Ly,rho,pi/2);

% Run the spin echo
report(spin_system,'eseem: running evolution to half-time point...');
rho_stack=evolution(spin_system,L,[],rho,parameters.timestep/2,(parameters.npoints-1),'trajectory',screen);

report(spin_system,'eseem: applying the refocusing pulse...');
rho_stack=step(spin_system,Ly,rho_stack,pi);

report(spin_system,'eseem: running evolution to refocus...');
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep/2,[],'refocus',coil);

% Detect
fid=transpose(coil'*rho_stack);

end
