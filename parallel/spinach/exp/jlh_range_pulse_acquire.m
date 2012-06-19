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

function fid=jlh_range_pulse_acquire(spin_system,parameters,L,rho)

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
            whosQsub = whos('Q[[1]]');
            whosIso = whos('Iso');
            msg = sprintf('jlh_pulse_acquire: dim(Iso) = %d x %d = %d ~ %d bytes ; dim(Q) = %d x %d = %d ~ %d bytes', whosIso.size(1), whosIso.size(2), whosIso.size(1) * whosIso.size(2), whosIso.bytes, whosQ.size(1), whosQ.size(2), whosQ.size(1) * whosQ.size(2), whosQ.bytes);
            report(spin_system,msg);
            
            msg = strcat('jlh_pulse_acquire: class(Iso) = ', whosIso.class, ' ; class(Q) = ', whosQ.class);
            report(spin_system,msg);
            
            %msg = strcat(    sprintf('jlh_pulse_acquire: dim(Q[]) = %d x %d = %d ~ %d bytes ; ', whosQsub.size(1), whosQsub.size(2), whosQsub.size(1) * whosQsub.size(2), whosQsub.bytes),' class(Q[]) = ', whosQsub.class );
            %report(spin_system,msg);
            
            
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
            %fid=zeros(parameters.npoints,1);
            
            % Loop over orientations
%jlh -------------------------------------------------
            orientations_per_node = floor(grid_size / parameters.nodes);
            remaining_orientations = mod(grid_size,parameters.nodes);
            %report(spin_system, strcat('jlh_pulce_acquire: ', grid_size, ' orientation / ', nodes, ' nodes = ',orientations_per_node, ' orientations per node; ', remaining_orientations, ' remaining orientations.'));
            report(spin_system, sprintf('jlh_pulse_acquire: %d orientation / %d nodes = %d orientations per node; %d remaining orientations.',grid_size,parameters.nodes,orientations_per_node,remaining_orientations));

            right_now = datestr(now,'yyyymmdd_HHMMSS');
            dir_name = [  func2str(@jlh_range_pulse_acquire) '_' right_now ];
            mkdir(dir_name);
            
            input_file_name = [ dir_name '/' func2str(@jlh_range_pulse_acquire) '_' right_now  '.mat'];
            report(spin_system, strcat('jlh_pulse_acquire: Saving spin_system, Iso, Q, Ly, rho, phi, theta, weight, coil, parameters into "', input_file_name, '"...'));
            save(input_file_name, 'spin_system','Iso', 'Q', 'Ly', 'rho', 'phi', 'theta', 'weight', 'coil','parameters');
            
            job_name = cell(parameters.nodes,1);
            job_identifier = cell(parameters.nodes,1);
            jobfile_name = cell(parameters.nodes,1);
            logfile_name = cell(parameters.nodes,1);
            output_file_name = cell(parameters.nodes,1);

            lower_boundary_offset =  [ zeros(1,parameters.nodes-remaining_orientations+1) (1:(remaining_orientations-1)) ];
            upper_boundary_offset = [ zeros(1,parameters.nodes-remaining_orientations) (1:remaining_orientations) ];
            lower_boundary = ( ((1:parameters.nodes)-1)*orientations_per_node+1 ) + lower_boundary_offset;
            upper_boundary = (1:parameters.nodes)*(orientations_per_node) + upper_boundary_offset;
            
            
            
            
            for i=1:parameters.nodes
                job_name(i) = cellstr( sprintf( strcat('jlh_', right_now, '_%d'), i) );
                jobfile_name(i) = cellstr( [ dir_name '/' job_name{i} '-jobfile'] );
                logfile_name(i) = cellstr( [ dir_name '/' job_name{i} '.log'] );
                output_file_name(i) = cellstr( [ dir_name '/' job_name{i} '.mat']);
                %lower_boundary = (i-1)*orientations_per_node+1;
                %upper_boundary = i*(orientations_per_node);
                %job_name = strcat('jlh_', right_now, '_', i);
                %output_file_name = strcat(job_name,'.mat');
                
                report(spin_system, sprintf('jlh_pulse_acquire: Preparing jlh_outsourced_pulse_acquire instance %d of %d...', i, parameters.nodes));
                report(spin_system, sprintf('jlh_pulse_acquire: Orientation %d to %d, in total %d orientations',lower_boundary(i),upper_boundary(i),upper_boundary(i)-lower_boundary(i)+1));

            %jlh - writing a jobfile for every set to hand to qsub: 
                command = [sprintf('/net/opt/bin/matlab -nodesktop -nosplash -r "jlh_outsourced_pulse_acquire(%d,%d',lower_boundary(i),upper_boundary(i)) ',''' input_file_name ''',''' output_file_name{i} '''); quit;" &>' logfile_name{i}];
        
                fileID = fopen(jobfile_name{i},'w');
                fprintf(fileID, '#!/bin/bash \n');
                fprintf(fileID, ['#PBS -N ' job_name{i} ' \n']);
                fprintf(fileID, '#PBS -q batch \n');
                fprintf(fileID, '#PBS -l nodes=1:ppn=%d \n', parameters.ppn);
                fprintf(fileID, '#PBS -cdl walltime=00:10:00 \n');
                fprintf(fileID, '#PBS -m bea -M jotelha@zedat.fu-berlin.de \n\n');
                
                fprintf(fileID, 'cd $PBS_O_WORKDIR\n\n');

                fprintf(fileID, [ command '\n']);
                
                fclose(fileID);
                
            %output the content of the jobfile
                jobfile_content = fileread(jobfile_name{i});
                report(spin_system, ['jlh_pulse_acquire: Creating jobfile "' jobfile_name{i} '": ' jobfile_content ]);

                %qsub_call = [sprintf('qsub -q batch -l nodes=1:ppn=%d',ppn),' -l walltime=00:10:00 -m bae -M jotelha@zedat.fu-berlin.de -N ' job_name{i} ' ' command];
                qsub_call = ['qsub ' jobfile_name{i}];
                
                report(spin_system, strcat('jlh_pulse_acquire: Calling "', qsub_call, '"...'));
                
                [status,stdout] = unix(qsub_call);
                if status~=0
                    report(spin_system, strcat('jlh_pulse_acquire: Error executing qsub: ', stdout ) );
                    return;
                end
                job_identifier(i) = cellstr(stdout);
                report(spin_system, strcat('jlh_pulce_acquire: Job "', job_identifier{i} , '" queued.'));
            end
            
            
%           last_set_starting_point = parameters.nodes*orientations_per_node+1;
%            report(spin_system, sprintf('jlh_pulce_acquire: Calculating remaining orientations %d to %d...', last_set_starting_point, grid_size ));    
            
            
%             parforTic = tic;
%             parfor n=last_set_starting_point:grid_size %jlh only certain range
%                 %jlh - benchmarking
%                 singleOrientationTic = tic;
%                 msg1 = sprintf('jlh_pulce_acquire: Orientation %d', n);
%                    
%                 % Get the current orientation
%                 % jlh - the only line differing from for to for
%                 L=Iso+orientation(Q,[phi(n) theta(n) 0]);
%                 
%                 %jlh - report the size of the Louvillian
%                 sizeL = size(L);
%                 s1 = sizeL(1);
%                 s2 = sizeL(2);
%                 stot = s1*s2;
%                 mbytes = double(stot*8)/(1024*1024);
%                 msg2 =  [ sprintf('jlh_pulse_acquire: dim(L) = %d x %d = %d ~ %.2f MB', s1,s2,stot,mbytes) ' ; class(L) = ' class(L) ];
%                 
%                 % Apply the pulse
%                 rho_pulsed=step(spin_system,Ly,rho,pi/2);
%             
%                 % Run the simulation
%                 fid=fid+weight(n)*evolution(spin_system,L,coil,rho_pulsed,timestep,parameters.npoints-1,'observable'); %#ok<PFBNS>
%                 
%                 singleOrientationTime = toc(singleOrientationTic);
%                 msg3 = sprintf('jlh_pulse_acquire: orientation computation time: %f s', singleOrientationTime);
%                 report(spin_system, msg1);
%                 report(spin_system, msg2);
%                 report(spin_system, msg3);
%             end 
%             parforTime = toc(parforTic);
%             
%             timeSpentMsg = sprintf('jlh_pulse_acquire: parfor computation time: %f s', parforTime);
%             report(spin_system, timeSpentMsg);
            
            stopped = false(parameters.nodes,1);
            no_of_stopped = 0;
            allStopped = false;
            
            pause on;
            while ~allStopped
                allStopped = true;
                for i=1:parameters.nodes
                    if ~stopped(i)
                        qstat_call = ['qstat ' job_identifier{i}];

                        report(spin_system, ['jlh_pulse_acquire: Calling "' qstat_call '"...']);

                        [status,stdout] = unix(qstat_call);
                        report(spin_system, sprintf('jlh_pulse_acquire: qstat output length: %d', length(stdout) ));
                        report(spin_system, ['jlh_pulse_acquire: qstat output:' stdout]);

                        if status==153
                            stopped(i) = true;
                            report(spin_system, ['jlh_pulse_acquire: Job "' job_identifier{i} '" stopped.']);
                            no_of_stopped = no_of_stopped + 1;
                        elseif status==0
                            report(spin_system, ['jlh_pulse_acquire: Job "' job_identifier{i} '" running.']);
                        else
                            report(spin_system, [sprintf('jlh_pulse_acquire: Error code (%d) executing qstat: ',status) stdout] );
                            return;
                        end
                    end
                    allStopped = allStopped && stopped(i);
                end
                report(spin_system, sprintf('jlh_pulse_acquire: %d Jobs stopped. Waiting 10s...', no_of_stopped));
                pause(10);
            end
            pause off;
            
            banner(spin_system,'outsourced_pulse_acquire_banner');
            fid=zeros(parameters.npoints,1);
            
            for i=1:parameters.nodes
                if exist(output_file_name{i}, 'file')
                    tmp = load( output_file_name{i}, 'fid');
                    fid = fid + tmp.fid;
                else
                   report(spin_system, [sprintf('jlh_pulse_acquire: Error: jlh_outsourced_pulse_acquire instance %d did not write output file "', i) output_file_name{i} '"!']);
                end
       %jlh - print jlh_outsourced_pulse_acquire output on stdout
                logfile_content = fileread(logfile_name{i});
                fprintf(logfile_content);
            end
    end
end