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

function fid=jlh_1parfor_pulse_acquire(spin_system,parameters,L,rho)

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
    
    % Run the simulation
    fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
            
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
            
            % Run the simulation
            fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
            
        case 'crystal'
            
            % Get the anisotropic Liouvillian
            report(spin_system,'pulse_acquire: building anisotropic Liouvillian...');
            [Iso,Q]=h_superop(spin_system); L=Iso+orientation(Q,[0 0 0]);
            
            %jlh - report data to be copied to the machines
            whosQ = whos('Q');
            whosIso = whos('Iso');
            msg = sprintf('jlh_pulse_acquire: dim(Iso) = %d x %d = %d ~ %.2f ; dim(Q) = %d x %d = %d ~ %.2f', whosIso.size(1), whosIso.size(2), whosIso.size(1) * whosIso.size(2), double(whosIso.bytes) / (1024*1024), whosQ.size(1), whosQ.size(2), whosQ.size(1) * whosQ.size(2), double(whosQ.bytes) / (1024*1024));
            report(spin_system,msg);
            
            % Apply the offset
            if isfield(parameters,'offset')
                report(spin_system,'pulse_acquire: applying the offset...');
                L=L-offset(spin_system,parameters.spins,parameters.offset);
            end
            
            % Get the detection state
            coil=state(spin_system,'L+',parameters.spins);
                        
            % Apply the pulse
            rho=step(spin_system,Ly,rho,pi/2);
            
            % Run the simulation
            fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
        
        case 'powder'
            
            % Get the anisotropic Liouvillian
            report(spin_system,'pulse_acquire: building anisotropic Liouvillian...');
            [Iso,Q]=h_superop(spin_system);
            
             %jlh - report data to be copied to the machines
            whosQ = whos('Q');
            whosIso = whos('Iso');
            msg = sprintf('jlh_pulse_acquire: dim(Iso) = %d x %d = %d ~ %.2f MB; dim(Q) = %d x %d = %d ~ %.2f MB', whosIso.size(1), whosIso.size(2), whosIso.size(1) * whosIso.size(2), double(whosIso.bytes) / (1024*1024), whosQ.size(1), whosQ.size(2), whosQ.size(1) * whosQ.size(2), double(whosQ.bytes) / (1024*1024));
            report(spin_system,msg);
            
            % Apply the offset
            if isfield(parameters,'offset')
                report(spin_system,'pulse_acquire: applying the offset...');
                Iso=Iso-offset(spin_system,parameters.spins,parameters.offset);
            end
    
            % Get the spherical averaging grid
            grid=load([spin_system.sys.root_dir spin_system.sys.slash 'exp' ...
                                                spin_system.sys.slash 'grids' ...
                                                spin_system.sys.slash 'lebedev_rank_' num2str(spin_system.tols.grid_rank) '.dat'],'ASCII');
            grid_size=size(grid,1); phi=pi*grid(:,1)/180; theta=pi*grid(:,2)/180; weight=grid(:,3);
            
            % Get the detection state
            coil=state(spin_system,'L+',parameters.spins);
            
            % Preallocate the answer
            fid=zeros(parameters.npoints,1);
            
            % Loop over orientations
            parforTic = tic;
            parfor n=1:1 %jlh only 1 repetitions
                %jlh - benchmarking
                singleOrientationTic = tic;
                msg1 = sprintf('jlh_pulce_acquire: Orientation %d', n);
                %report(spin_system, msg);
                
                % Get the current orientation
                L=Iso+orientation(Q,[phi(n) theta(n) 0]);
                
                %jlh - report the size of the Louvillian
                sizeL = size(L);
                s1 = sizeL(1);
                s2 = sizeL(2);
                stot = s1*s2;
                mbytes = double(stot*8)/(1024*1024);
                msg2 =  [ sprintf('jlh_pulse_acquire: dim(L) = %d x %d = %d ~ %.2f MB', s1,s2,stot,mbytes) ' ; class(L) = ' class(L) ];
                %msg = sprintf('jlh_pulse_acquire: dim(L) = %d x %d ; size(L) = %d', infoL.size(1), infoL.size(2), infoL.bytes);
                %report(spin_system,msg);
                
                % Apply the pulse
                rho_pulsed=step(spin_system,Ly,rho,pi/2);
            
                % Run the simulation
                fid=fid+weight(n)*evolution(spin_system,L,coil,rho_pulsed,timestep,parameters.npoints-1,'observable'); %#ok<PFBNS>
                
                singleOrientationTime = toc(singleOrientationTic);
                msg3 = sprintf('jlh_pulse_acquire: orientation computation time: %f s', singleOrientationTime);
                report(spin_system, msg1);
                report(spin_system, msg2);
                report(spin_system, msg3);
            end
            parforTime = toc(parforTic);
            timeSpentMsg = sprintf('parfor computation time: %f s', parforTime);
            report(spin_system, timeSpentMsg);
           
            
    end
end

end

