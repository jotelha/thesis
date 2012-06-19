% The standard 90-acqure pulse sequence, modified by jlh.
% Splits and parallelizes the propagation into equally sized packages of
% orientations.
%
% The following parameters are currently accepted:
%
%       parameters.spins      - '1H', '13C', etc
%       parameters.sweep       - the width of the spectral window (Hz)
%       parameters.npoints     - number time steps in the simulation
%       parameters.zerofill    - number of points to zerofill to
%       parameters.offset      - transmitter offset (Hz)
%       parameters.axis_units
%   Modification by jlh:
%       parameters.title       - the title of the simulation
%       parameters.nodes       - the no. of cluster nodes to acquire
%       parameters.ppn         - processors per node
%
% If a Liouvillian L is supplied, it is used for propagation, otherwise it
% would be generated here. If a state vector rho is supplied, it is used 
% as a starting point, otherwise thermal equilibrium is used.
%       
% ilya.kuprov@oerc.ox.ac.uk

function jlh_master_pulse_acquire(spin_system,parameters,L,rho)

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

% jlh - only allow the 'powder' simulation type

% Assemble the Liouvillian
% if exist('L','var')&&(~isempty(L))
%     
%     % Inform the user
%     report(spin_system,'pulse_acquire: using the Liouvillian as supplied.');
%     
%     % Apply the offset
%     if isfield(parameters,'offset')
%         report(spin_system,'pulse_acquire: applying the offset...');
%         L=L-offset(spin_system,parameters.spins,parameters.offset);
%     end
%     
%     % Get the detection state
%     coil=state(spin_system,'L+',parameters.spins);
%         
%     % Apply the pulse
%     rho=step(spin_system,Ly,rho,pi/2);
%     
%     % Run the simulation
%     fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
%             
% else
    
    switch spin_system.inter.regime
    
%         case 'liquid'
%             
%             % Get the isotropic Liouvillian
%             report(spin_system,'pulse_acquire: building isotropic Liouvillian...');
%             L=h_superop(spin_system)+1i*r_superop(spin_system)+1i*k_superop(spin_system);
%             
%             % Apply the offset
%             if isfield(parameters,'offset')
%                 report(spin_system,'pulse_acquire: applying the offset...');
%                 L=L-offset(spin_system,parameters.spins,parameters.offset);
%             end
%             
%             % Get the detection state
%             coil=state(spin_system,'L+',parameters.spins);
%                         
%             % Apply the pulse
%             rho=step(spin_system,Ly,rho,pi/2);
%             
%             % Run the simulation
%             fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
%             
%         case 'crystal'
%             
%             % Get the anisotropic Liouvillian
%             report(spin_system,'pulse_acquire: building anisotropic Liouvillian...');
%             [Iso,Q]=h_superop(spin_system); L=Iso+orientation(Q,[0 0 0]);
%             
%             % Apply the offset
%             if isfield(parameters,'offset')
%                 report(spin_system,'pulse_acquire: applying the offset...');
%                 L=L-offset(spin_system,parameters.spins,parameters.offset);
%             end
%             
%             % Get the detection state
%             coil=state(spin_system,'L+',parameters.spins);
%                         
%             % Apply the pulse
%             rho=step(spin_system,Ly,rho,pi/2);
%             
%             % Run the simulation
%             fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');
        
        case 'powder'
            
            % Get the anisotropic Liouvillian
            report(spin_system,'pulse_acquire: building anisotropic Liouvillian...');
            [Iso,Q]=h_superop(spin_system);
            
            %jlh - report size of isotropic Louvillian
            whosQ = whos('Q');
            whosIso = whos('Iso');
            msg = sprintf('jlh_pulse_acquire: dim(Iso) = %d x %d = %d ~ %d bytes ; dim(Q) = %d x %d = %d ~ %d bytes', whosIso.size(1), whosIso.size(2), whosIso.size(1) * whosIso.size(2), whosIso.bytes, whosQ.size(1), whosQ.size(2), whosQ.size(1) * whosQ.size(2), whosQ.bytes);
            report(spin_system,msg);
            
            msg = strcat('jlh_pulse_acquire: class(Iso) = ', whosIso.class, ' ; class(Q) = ', whosQ.class);
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
            
           
%jlh --------------------------------------------------
%jlh Arranging orientation packages and queuing threads
%jlh --------------------------------------------------

            % calculating the integer part of orientations per node package
            orientations_per_node = floor(grid_size / parameters.nodes);
            % calculating the remaining orientations to be distributed
            remaining_orientations = mod(grid_size,parameters.nodes);
            report(spin_system, sprintf('jlh_pulse_acquire: %d orientation / %d nodes = %d orientations per node; %d remaining orientations.',grid_size,parameters.nodes,orientations_per_node,remaining_orientations));

            % preparing filenames:
            
            % files named with the title of the simulation
            % summary_jobfile_name  -- Jobfile for summary process jlh_sum_results() to be run eventually collecting all data
            % summary_logfile_name  -- Logfile of every node package propagation redirected to here
            % spectrum_file_name    -- .dat fiel for the resulting spectrum
            % input_file_name       -- .mat variables dump to be handed to the node package processes
            
            % files named with the title of the simulation and integer index
            % job_name              -- names of the prepared Jobs to be written to the jobfiles
            % job_identifier        -- the job ID returned by sbatch when qeueing
            % jobfile_name          -- Jobfiles for jlh_outsourced_pulse_acquire() node package processes
            % logfile_name          -- Logfiles of every node package process to be collected by summary process
            base_name = parameters.title;
            %creating an own directory to collect all the data of the
            %current simulation:
            dir_name = base_name;
            mkdir(dir_name);
            
            summary_jobfile_name = [ dir_name '/'  base_name '-jobfile'];
            summary_logfile_name = [ dir_name '/'  base_name '.log'];
            spectrum_file_name = [ dir_name '/' base_name '.dat'];
             
            input_file_name = [ dir_name '/' base_name '.mat'];
            
            job_name = cell(parameters.nodes,1);
            job_identifier = cell(parameters.nodes,1);
            jobfile_name = cell(parameters.nodes,1);
            logfile_name = cell(parameters.nodes,1);
            output_file_name = cell(parameters.nodes,1);
            for i=1:parameters.nodes
                job_name(i) = cellstr( sprintf( strcat(base_name,'_%d'), i) );
                jobfile_name(i) = cellstr( [ dir_name '/' job_name{i} '-jobfile'] );
                logfile_name(i) = cellstr( [ dir_name '/' job_name{i} '.log'] );
                output_file_name(i) = cellstr( [ dir_name '/' job_name{i} '.mat']);
            end

            %Distributing the orientations on the nodes and adding
            %remaining orientations to some packages:
            if remaining_orientations == 0 
                lower_boundary_offset = zeros(1,parameters.nodes);
                upper_boundary_offset = zeros(1,parameters.nodes);
            else
                lower_boundary_offset =  [ zeros(1,parameters.nodes-remaining_orientations+1) (1:(remaining_orientations-1)) ];
                upper_boundary_offset = [ zeros(1,parameters.nodes-remaining_orientations) (1:remaining_orientations) ];
            end
            lower_boundary = ( ((1:parameters.nodes)-1)*orientations_per_node+1 ) + lower_boundary_offset;
            upper_boundary = (1:parameters.nodes)*(orientations_per_node) + upper_boundary_offset;
            
            report(spin_system, strcat('jlh_pulse_acquire: Saving spin_system, Iso, Q, Ly, rho, phi, theta, weight, coil, parameters into "', input_file_name, '"...'));
            %Writing variables dump to .mat-file
            save(input_file_name, 'spin_system','Iso', 'Q', 'Ly', 'rho', 'phi', 'theta', 'weight', 'coil','parameters','output_file_name', 'logfile_name', 'spectrum_file_name');
            
         
            %writing jobfiles for every node package process, evoking a
            %seperate matlab process calling jlh_outsourced_pulse_acquire()
            %each
            for i=1:parameters.nodes    
                report(spin_system, sprintf('jlh_pulse_acquire: Preparing jlh_outsourced_pulse_acquire instance %d of %d...', i, parameters.nodes));
                report(spin_system, sprintf('jlh_pulse_acquire: Orientation %d to %d, in total %d orientations',lower_boundary(i),upper_boundary(i),upper_boundary(i)-lower_boundary(i)+1));

                %the matlab command to be run
                command = [sprintf('matlab -nodesktop -nosplash -r "jlh_outsourced_pulse_acquire(%d,%d',lower_boundary(i),upper_boundary(i)) ',''' input_file_name ''',''' output_file_name{i} '''); quit;" &>' logfile_name{i}];
        
                fileID = fopen(jobfile_name{i},'w');
                %setting SBATCH options
                fprintf(fileID, '#!/bin/bash \n');
                fprintf(fileID, ['#SBATCH --job-name=' job_name{i} ' \n']);
                %fprintf(fileID, '#PBS -q batch \n');
                fprintf(fileID, '#SBATCH --nodes=1 \n');
                fprintf(fileID, '#SBATCH --ntasks-per-node=%d \n', parameters.ppn);
                %fprintf(fileID, '#SBATCH --time=20:00:00 \n');
                fprintf(fileID, '#SBATCH --mail-type=ALL\n');
                fprintf(fileID, '#SBATCH --mail-user=jotelha@zedat.fu-berlin.de \n');
                %redirecting torque logfiles to simulation directory
                fprintf(fileID, ['#SBATCH --error=' dir_name '/' job_name{i} '.e \n']);
                fprintf(fileID, ['#SBATCH --output=' dir_name '/' job_name{i} '.o \n\n']);
                
                fprintf(fileID, 'cd $SUBMITDIR\n\n');
                fprintf(fileID, 'module load matlab/R2011b\n\n');

                fprintf(fileID, [ command '\n']);
                
                fclose(fileID);
                
                %output the content of the jobfile
                jobfile_content = fileread(jobfile_name{i});
                report(spin_system, ['jlh_pulse_acquire: Creating jobfile "' jobfile_name{i} '": ' jobfile_content ]);
                
                %sbatch command to be run by unix shell
                sbatch_call = ['sbatch ' jobfile_name{i}];
                
                report(spin_system, strcat('jlh_pulse_acquire: Calling "', sbatch_call, '"...'));
                
                %calling sbatch via unix interface
                [status,stdout] = unix(sbatch_call);
                
                %checking for sbatch return status 0 => job queued
                %successfully
                if status~=0
                    report(spin_system, strcat('jlh_pulse_acquire: Error executing sbatch: ', stdout ) );
                    return;
                end
                %job identifier extracted from sbatch stdout
                job_identifier(i) = cellstr( num2str( sscanf(stdout, 'Submitted batch job %d') ) );
                report(spin_system, strcat('jlh_pulce_acquire: Job "', job_identifier{i} , '" queued.'));
            end
            
            
            %writing a jobfile for summary process depending on the
            %successfull completion of all other node package processes,
            %structure similar to above
           
            command = ['matlab -nodesktop -nosplash -r "jlh_sum_results(''' input_file_name '''); quit;" &> ' summary_logfile_name];

            fileID = fopen(summary_jobfile_name,'w');
            %fprintf(fileID, '#!/bin/bash \n');
            %fprintf(fileID, ['#PBS -N ' base_name ' \n']);
            %fprintf(fileID, '#PBS -q batch \n');
            %fprintf(fileID, '#PBS -l nodes=1:ppn=1\n');
            %fprintf(fileID, '#PBS -l walltime=01:00:00 \n');
            %fprintf(fileID, '#PBS -m bea -M jotelha@zedat.fu-berlin.de \n');
            %fprintf(fileID, ['#PBS -e ' base_name '\n']);
            %fprintf(fileID, ['#PBS -o ' base_name '\n']);
            
            fprintf(fileID, '#!/bin/bash \n');
            fprintf(fileID, ['#SBATCH --job-name=' base_name  '\n']);
            %fprintf(fileID, '#PBS -q batch \n');
            fprintf(fileID, '#SBATCH --nodes=1 \n');
            fprintf(fileID, '#SBATCH --ntasks-per-node=1 \n');
            %fprintf(fileID, '#SBATCH --time=1:00:00 \n');
            fprintf(fileID, '#SBATCH --mail-type=ALL\n');
            fprintf(fileID, '#SBATCH --mail-user=jotelha@zedat.fu-berlin.de \n');
            %redirecting torque logfiles to simulation directory
            fprintf(fileID, ['#SBATCH --error=' dir_name '/' base_name '.e \n']);
            fprintf(fileID, ['#SBATCH --output=' dir_name '/' base_name '.o \n']);
                
            %formatting job dependencies list 'afterok:job1:job2:...:jobN'
            tmpstr = strcat(':',job_identifier);
            dependency_list = [tmpstr{:}];
            
            fprintf(fileID, ['#SBATCH --dependency=afterok' dependency_list '\n\n']);
            
            fprintf(fileID, 'cd $SUBMITDIR\n\n');
            fprintf(fileID, 'module load matlab/R2011b\n\n');

            fprintf(fileID, [ command '\n']);

            fclose(fileID);
            
            summary_jobfile_content = fileread(summary_jobfile_name);
            report(spin_system, ['jlh_pulse_acquire: Creating jobfile "' summary_jobfile_name '": ' summary_jobfile_content ]);

            sbatch_call = ['sbatch ' summary_jobfile_name];

            report(spin_system, strcat('jlh_pulse_acquire: Calling "', sbatch_call, '"...'));

            %queing summary process
            [status,stdout] = unix(sbatch_call);
            if status~=0
                report(spin_system, strcat('jlh_pulse_acquire: Error executing sbatch: ', stdout ) );
                return;
            end
            
            report(spin_system, strcat('jlh_pulce_acquire: ', stdout));          
    end
end