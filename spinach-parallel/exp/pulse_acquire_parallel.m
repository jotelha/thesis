% The standard 90-acqure pulse sequence.
%
% The following parameters are currently accepted:
%
%       parameters.spins      - '1H', '13C', etc
%       parameters.sweep       - the width of the spectral window (Hz)
%       parameters.npoints     - number time steps in the simulation
%       parameters.zerofill    - number of points to zerofill to
%       parameters.offset      - transmitter offset (Hz)
%       parameters.axis_units
%
% If a Liouvillian L is supplied, it is used for propagation, otherwise it
% would be generated here. If a state vector rho is supplied, it is used 
% as a starting point, otherwise thermal equilibrium is used.
%       
% ilya.kuprov@oerc.ox.ac.uk

function fid=pulse_acquire_parallel(spin_system,parameters,L,rho)

% Show the banner
banner(spin_system,'sequence_banner');

% Print a report
report(spin_system,'pulse_acquire: computing pulse_acquire...');
sequence_report(spin_system,parameters);

% Compute the digitization parameters.
timestep=1/parameters.sweep;

% Generate the basic operators
Lp=operator(spin_system,'L+',parameters.spins);
Ly=(Lp-Lp')/2i;

% Set the secularity assumptions
spin_system=secularity(spin_system,'nmr');

% Get the initial state
if exist('rho','var')&&(~isempty(rho))
    
    % Inform the user
    report(spin_system,'pulse_acquire: using the initial state vector as supplied.');
    
else
    
    % Start from thermal equilibrium
    rho=equilibrium(spin_system);
    
end

% Assemble the Liouvillian
if exist('L','var')&&(~isempty(L))
    
    % Inform the user
    report(spin_system,'pulse_acquire: using the Liouvillian as supplied.');
    
    % Apply the offset
    if isfield(parameters,'offset')
        report(spin_system,'pulse_acquire: applying the offset...');
        L=L-offset(spin_system,parameters.spins,parameters.offset);
    end
    
    % Get the detection state
    coil=state(spin_system,'L+',parameters.spins);
        
    % Apply the pulse
    rho=step(spin_system,Ly,rho,pi/2);
            
else
    
    switch spin_system.inter.regime
    
        case 'liquid'
            
            % Get the isotropic Liouvillian
            report(spin_system,'pulse_acquire: building isotropic Liouvillian...');
            L=h_superop(spin_system)+1i*r_superop(spin_system)+1i*k_superop(spin_system);
            
            % Apply the offset
            if isfield(parameters,'offset')
                report(spin_system,'pulse_acquire: applying the offset...');
                L=L-offset(spin_system,parameters.spins,parameters.offset);
            end
            
            % Get the detection state
            coil=state(spin_system,'L+',parameters.spins);
                        
            % Apply the pulse
            rho=step(spin_system,Ly,rho,pi/2);
            
        case 'crystal'
            
            % Get the anisotropic Liouvillian
            report(spin_system,'pulse_acquire: building anisotropic Liouvillian...');
            [Iso,Q]=h_superop(spin_system); L=Iso+orientation(Q,[0 0 0]);
            
            % Apply the offset
            if isfield(parameters,'offset')
                report(spin_system,'pulse_acquire: applying the offset...');
                L=L-offset(spin_system,parameters.spins,parameters.offset);
            end
            
            % Get the detection state
            coil=state(spin_system,'L+',parameters.spins);
                        
            % Apply the pulse
            rho=step(spin_system,Ly,rho,pi/2);
        
        case 'powder'
            
            % Get the anisotropic Liouvillian
            report(spin_system,'pulse_acquire: building anisotropic Liouvillian...');
            [Iso,Q]=h_superop(spin_system);
            
            % Apply the offset
            if isfield(parameters,'offset')
                report(spin_system,'pulse_acquire: applying the offset...');
                Iso=Iso-offset(spin_system,parameters.spins,parameters.offset);
            end
    
            % Get the spherical averaging grid
            grid=load([spin_system.sys.root_dir spin_system.sys.slash 'exp' ...
		       spin_system.sys.slash 'grids' spin_system.sys.slash ...
		       'lebedev_rank_' num2str(spin_system.tols.grid_rank) ...
		       '.dat'],'ASCII');
	    
	    % Preallocate the answer
	    fid=zeros(parameters.npoints,1);
	    temp_fid=zeros(parameters.npoints,1);
	    %How many orientations assigned to one core??
	    %Tuned to 8GB RAM/core for 20 spins
	    gridstep=347;
	    
	    % Initial state
	    rho=kron(ones(gridstep,1),rho);
	    
	    parfor n=1:10
	      gridi=(n-1)*gridstep+1;
	      gridf=n*gridstep;
	      phi=pi*grid(gridi:gridf,1)/180;
	      theta=pi*grid(gridi:gridf,2)/180;
	      weight=grid(gridi:gridf,3);

	      % Get the orientation array
	      L_aniso=orientation(Q,[phi theta zeros(size(theta))]);
	      L=blkdiag(L_aniso{:})+kron(speye(gridstep),Iso);
	      L=clean_up(spin_system,L,spin_system.tols.liouv_zero);
	      
	      %Detection state, apply grid weights
	      coil=kron(ones(gridstep,1),state(spin_system,'L+',parameters.spins));
	      %coil=kron(weight,state(spin_system,'L+',parameters.spins));
	      
	      % Apply the pulse
	      rho_pulsed=step(spin_system,kron(diag(weight),Ly),rho,pi/2);
	      
	      % Run the simulation
	      temp_fid=evolution(spin_system,L,coil,rho_pulsed,timestep, ...
				 parameters.npoints-1,'observable');
	      fid=fid+temp_fid;
    end
end
end
