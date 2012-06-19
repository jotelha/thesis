% Mims and continuous-wave ENDOR pulse sequences.
%
% The following parameters are accepted:
%
% parameters.type               specifying 'mims' produces Mims ENDOR,
%                               specifying 'cw' produces continuous-
%                               wave ENDOR
%
% parameters.sweep              nuclear frequency sweep width, Hz
%
% parameters.npoints            number of fid points to be computed
%
% parameters.tau                stimulated echo time for Mims ENDOR,
%                               seconds
%
% parameters.electron_spin      set to 'E' to pulse all electrons in the
%                               system, specify a single number to pulse
%                               a specific electron
%
% parameters.nuclear_spins      nuclei to be detected
%
% parameters.orientation        three Euler angles (column vector, in
%                               radians) specifying molecular orientation
%                               relative to the input orientation if
%                               'crystal' was chosen in inter.regime
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function fid=endor(spin_system,parameters)

% Show the banner
banner(spin_system,'sequence_banner');

%% Build the simulation infrastructure

% Print a report
report(spin_system,'endor: computing ENDOR...');
sequence_report(spin_system,parameters);

% Set the sequence-specific parameters
switch parameters.type
    
    case 'mims'
        
        % Generate electron and nuclear pulse operators
        Ly=(operator(spin_system,'L+',parameters.electron_spin)-...
            operator(spin_system,'L-',parameters.electron_spin))/2i;
        Sy=(operator(spin_system,'L+',parameters.nuclear_spins)-...
            operator(spin_system,'L-',parameters.nuclear_spins))/2i;
        
        % Detect the electron
        coil=state(spin_system,'L+',parameters.electron_spin);
        
        % Set the initial state to thermal equilibrium
        rho=equilibrium(spin_system);
        
    case 'cw'
        
        % Generate the nuclear pulse operator
        Sy=(operator(spin_system,'L+',parameters.nuclear_spins)-...
            operator(spin_system,'L-',parameters.nuclear_spins))/2i;
        
        % Detect the nuclei
        coil=state(spin_system,'L+',parameters.nuclear_spins);
        
        % Set the initial state to nuclear Lz, with populations proportional
        % to the hyperfine coupling
        rho=zeros(spin_system.bas.nstates,1);
        for L=find(strcmp(spin_system.comp.isotopes,parameters.electron_spin))
            for S=find(strcmp(spin_system.comp.isotopes,parameters.nuclear_spins))
                amplitude=norm(spin_system.inter.coupling.matrix{L,S})+norm(spin_system.inter.coupling.matrix{S,L});
                rho=rho+amplitude*state(spin_system,'Lz',S);
            end
        end
        rho=rho/norm(rho);
        
    otherwise
        
        error('endor: unknown simulation type.');

end
        
%% Generate the Liouvillian

% Set the secularity assumptions
spin_system=secularity(spin_system,'endor');

% Update the user
report(spin_system,'endor: building the Liouvillian...');

% Decide the simulation regime
switch spin_system.inter.regime
    
    case 'liquid'
        
        % Get the Liouvillian
        L=h_superop(spin_system)+1i*r_superop(spin_system)+1i*k_superop(spin_system);      
        
    case 'crystal'
        
        % Get the Liouvillian
        [Isotropic,RotationalBasis]=h_superop(spin_system);
        Anisotropic=orientation(RotationalBasis,parameters.orientation);
        L=Isotropic+Anisotropic{1}+1i*r_superop(spin_system)+1i*k_superop(spin_system);    
       
    case 'powder'
        
        % Get the Liouvillian
        [Isotropic,RotationalBasis]=h_superop(spin_system);
        
        % Get the spherical averaging grid
        grid=load([spin_system.sys.root_dir spin_system.sys.slash 'exp' ...
                                            spin_system.sys.slash 'grids' ...
                                            spin_system.sys.slash 'lebedev_rank_' num2str(spin_system.tols.grid_rank) '.dat'],'ASCII');
        grid_size=size(grid,1); phi=pi*grid(:,1)/180; theta=pi*grid(:,2)/180; weight=grid(:,3);
        
        % Get the orientation array for the Liouvillian
        L_aniso=orientation(RotationalBasis,[phi theta zeros(size(theta))]);
        L=blkdiag(L_aniso{:})+kron(speye(grid_size),Isotropic+1i*r_superop(spin_system)+1i*k_superop(spin_system));
        
        % GGet the orientation arrays for the rest of the infrastructure
        rho=kron(ones(grid_size,1),rho); coil=kron(weight,coil); 
        Ly=kron(speye(grid_size),Ly); Sy=kron(speye(grid_size),Sy);
       
end

%% Run the simulation

switch parameters.type
    
    case 'mims'
        
        % Choose an appropriate stepsize and number of steps for the tau evolution
        [timestep,nsteps]=stepsize(L,parameters.tau);
        
        % Apply the initial pulses
        rho=step(spin_system,Ly,rho,pi/2);
        rho=evolution(spin_system,L,[],rho,timestep,nsteps,'final');
        rho=step(spin_system,Ly,rho,pi/2);
        
        % Apply pulses on nuclear spins
        rho=step(spin_system,Sy,rho,pi/2);
        rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep,parameters.npoints-1,'trajectory');
        rho_stack=step(spin_system,Sy,rho_stack,pi/2)-step(spin_system,-Sy,rho_stack,pi/2);
        
        % Apply pulse on electron spin to refocus
        rho_stack=step(spin_system,Ly,rho_stack,pi/2);
        rho_stack=evolution(spin_system,L,[],rho_stack,timestep,nsteps,'final',coil);
        
        % Detect
        fid=transpose(coil'*rho_stack);
        
    case 'cw'
        
        % Compute the digitization parameters.
        timestep=1/parameters.sweep;
        
        % Run the simulation
        rho=step(spin_system,Sy,rho,pi/2);
        fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
        
        % Symmetrize the frequencies
        fid=real(fid);

end

end
