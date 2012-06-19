% Validation module. Inspects user input to various functions and complains
% to the console if something is not right.
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk
% matthew.krzystyniak@oerc.ox.ac.uk

function validate(spin_system,parameter,check)

%% Validation tolerances

% Maximum imaginary component in a coupling that is meant to be real
max_imag=1e-3;

% Maximum non-Hermiticity in a matrix that is meant to be Hermitian
max_nonherm=1e-3;

%% Validation checks

switch check
    
    case 'create_runtime_input'              % Done
        
        % Move to the natural name
        sys=parameter;
        
        % Check the hush switch
        if isfield(sys,'output')
            if ~any(strcmp(sys.output,{'console','file','hush'}))
                    error('validate: the allowed values for sys.output are `console`, `file` and `hush`.');
            end
        end

        % Check the log file redirect
        if isfield(sys,'logfile')
            if ~ischar(sys.logfile)
                error('validate: log file name in sys.logfile must be a string.');
            end
        end
        
        % Check the disable switches
        if isfield(sys,'disable')
            for n=1:numel(sys.disable)
                if ~any(strcmp(sys.disable{n},{'zte','pt','symmetry','krylov','clean-up','dss','expv','caching'}))
                    error('validate: the allowed values for sys.disable are `zte`, `pt`, `symmetry`, `krylov`, `clean-up`, `dss`, `expv`, and `caching`.');
                end
            end
        end
        
    case 'create_comp_input'                 % Done
        
        % Move to the natural name
        sys=parameter;
        
        % Check the caching identifier
        if isfield(sys,'identifier')
            if ~ischar(sys.identifier)
                error('validate: the caching identifier in sys.identifier must be a string.');
            end
        end
        
        % Check the regime switch
        if isfield(sys,'regime')
            if ~any(strcmp(sys.regime,{'liquid','crystal','powder'}))
                error('validate: the allowed values for sys.regime are `liquid`, `crystal` and `powder`.');
            end
        end
        
        % Check the isotopes variable
        if ~isfield(sys,'isotopes')
            error('validate: sys.isotopes variable must be declared.');
        elseif isempty(sys.isotopes)
            error('validate: sys.isotopes cannot be empty.');
        elseif ~iscell(sys.isotopes)
            error('validate: sys.isotopes must be a cell array.');
        elseif ~all(cellfun(@ischar,sys.isotopes))
            error('validate: all elements of sys.isotopes cell array must be strings.');
        end
        
        % Check the labels variable
        if isfield(sys,'labels')
            if ~iscell(sys.labels)
                error('validate: sys.labels must be a cell array.');
            elseif isempty(sys.labels)
                error('validate: sys.labels cannot be empty.');
            elseif ~all(cellfun(@ischar,sys.labels))
                error('validate: all elements of sys.labels cell array must be strings.');
            elseif numel(sys.labels)~=numel(sys.isotopes)
                error('validate: the length of sys.labels must be the same as the length of sys.isotopes.');
            end
        end
        
    case 'create_zeeman_input'               % Done
        
        % Move to the natural names
        sys=parameter{1};
        inter=parameter{2};
        
        % Check the magnet variable
        if isfield(sys,'magnet')
            if (~isnumeric(sys.magnet))||(~isreal(sys.magnet))||(numel(sys.magnet)~=1)
                error('validate: sys.magnet must be a real number.');
            end
        end
        
        % Check the inter.zeeman variable
        if isfield(inter,'zeeman')
            
            % Check eigenvalues / Euler angles specification
            if isfield(inter.zeeman,'eigs')
                
                % Check the type
                if ~iscell(inter.zeeman.eigs)
                    error('validate: inter.zeeman.eigs must be a cell array of 1x3 vectors.');
                end
                
                % Check the dimensions
                if numel(inter.zeeman.eigs)~=spin_system.comp.nspins
                    error('validate: the number of elements in inter.zeeman.eigs must match the number of spins.');
                end
                
                % Make sure eulers exist
                if ~isfield(inter.zeeman,'euler')
                    error('validate: inter.zeeman.euler variable must be set together with inter.zeeman.eigs.');
                end
                
                % Make sure eulers are cells
                if ~iscell(inter.zeeman.euler)
                    error('validate: inter.zeeman.euler must be a cell array of 1x3 vectors.');
                end
                
                % Make sure eulers have the correct length
                if ~all(size(inter.zeeman.eigs)==size(inter.zeeman.euler))
                    error('validate: inter.zeeman.eigs and inter.zeeman.euler variables must have the same dimension.');
                end
                
                % Make sure all non-empty elements are 3-vectors
                for n=1:spin_system.comp.nspins
                
                    % For eigenvalues
                    if (~all(size(inter.zeeman.eigs{n})==[1 3]))&&(~all(size(inter.zeeman.eigs{n})==[0 0]))||(~isnumeric(inter.zeeman.eigs{n}))
                        error('validate: non-empty elements of inter.zeeman.eigs must be 1x3 vectors.');
                    end
                    
                    % For Euler angles
                    if (~all(size(inter.zeeman.euler{n})==[1 3]))&&(~all(size(inter.zeeman.euler{n})==[0 0]))||(~isnumeric(inter.zeeman.euler{n}))
                        error('validate: non-empty elements of inter.zeeman.euler must be 1x3 vectors.');
                    end
                    
                end
                
            end
            
            % Check the matrix specification
            if isfield(inter.zeeman,'matrix')
                
                % Check the type
                if ~iscell(inter.zeeman.matrix)
                    error('validate: inter.zeeman.matrix must be a cell array of 3x3 matrices.');
                end
                
                % Check the length 
                if numel(inter.zeeman.matrix)~=spin_system.comp.nspins
                    error('validate: the number of elements in the inter.zeeman.matrix array should match the number of spins.');
                end
                
                % Make sure all non-empty elements are 3x3 matrices 
                for n=1:spin_system.comp.nspins
                    if ((~all(size(inter.zeeman.matrix{n})==[3 3]))&&(~all(size(inter.zeeman.matrix{n})==[0 0])))||(~isnumeric(inter.zeeman.matrix{n}))
                        error('validate: non-empty elements of inter.zeeman.matrix must be 3x3 matrices.');
                    end
                end
                
            end
            
            % Check the scalars
            if isfield(inter.zeeman,'scalar')
                
                 % Check the type
                if ~iscell(inter.zeeman.scalar)
                    error('validate: inter.zeeman.scalar must be a cell array of numbers.');
                end
                
                % Check the length
                if numel(inter.zeeman.scalar)~=spin_system.comp.nspins
                    error('validate: the number of elements in the inter.zeeman.scalar array should match the number of spins.');
                end
                
                % Make sure all non-empty elements are numbers 
                for n=1:spin_system.comp.nspins
                    if (((numel(inter.zeeman.scalar{n})~=1)&&(numel(inter.zeeman.scalar{n})~=0)))||(~isnumeric(inter.zeeman.scalar{n}))
                        error('validate: non-empty elements of inter.zeeman.scalar must be numbers.');
                    end
                end
                
            end
            
        end
        
    case 'create_coupling_input'             % Done
        
        % Move to the natural name
        inter=parameter;
        
        % Check the coordinates
        if isfield(inter,'coordinates')
            
            % Check the type
            if ~iscell(inter.coordinates)
                error('validate: inter.coordinates must be a cell array of 1x3 vectors.');
            end
            
            % Check the size
            if numel(inter.coordinates)~=spin_system.comp.nspins
                error('validate: the number of elements in inter.coordinates must match the number of spins.')
            end
            
            % Check the contents
            for n=1:spin_system.comp.nspins
                
                % Make sure we have 3-vectors
                if ((~all(size(inter.coordinates{n})==[1 3]))&&(~all(size(inter.coordinates{n})==[0 0])))||(~isnumeric(inter.coordinates{n}))
                    error('validate: non-empty elements of inter.coordinates must be 1x3 vectors.');
                end
                
                % Make sure there are no collisions
                for k=1:spin_system.comp.nspins
                    if (~isempty(inter.coordinates{n}))&&(~isempty(inter.coordinates{k}))&&(n~=k)
                        if norm(inter.coordinates{n}-inter.coordinates{k})<0.5
                            error('validate: unphysical proximity detected in inter.coordinates.');
                        end
                    end
                end
                
            end
            
        end
        
        % Check the order matrix
        if isfield(inter,'order_matrix')
            
            % Make sure we have a 3x3 matrix
            if (~isnumeric(inter.order_matrix))||(~all(size(inter.order_matrix)==[3 3]))||(~isreal(inter.order_matrix))
                error('validate: inter.order_matrix must be a 3x3 matrix of real numbers.');
            end
            
            % Make sure the matrix is traceless
            if trace(inter.order_matrix)>1e-10
                error('validate: inter.order_matrix must be traceless.');
            end
            
            % Make sure the matrix is symmetric
            if norm(inter.order_matrix-inter.order_matrix')>1e-10
                error('validate: inter.order_matrix must be symmetric.');
            end
            
            % Prohibit Redfield relaxation theory
            if isfield(inter,'relaxation')&&strcmp(inter.relaxation,'redfield')
                error('validate: Redfield theory is not available for liquid crystals.');
            end
            
        end
        
        % Check the couplings
        if isfield(inter,'coupling')
            
            % Check eigenvalues / Euler angles specification
            if isfield(inter.coupling,'eigs')
                
                % Check the type
                if ~iscell(inter.coupling.eigs)
                    error('validate: inter.coupling.eigs must be a cell array of 1x3 vectors.');
                end
                
                % Check the dimensions
                if ~all(size(inter.coupling.eigs)==[spin_system.comp.nspins spin_system.comp.nspins])
                    error('validate: both dimensions of inter.coupling.eigs must match the number of spins.');
                end
                
                % Make sure eulers exist
                if ~isfield(inter.coupling,'euler')
                    error('validate: inter.coupling.euler array must be set together with inter.coupling.eigs.');
                end
                
                % Make sure eulers are cells
                if ~iscell(inter.coupling.euler)
                    error('validate: inter.coupling.euler must be a cell array of 1x3 vectors.');
                end
                
                % Make sure eulers have the correct length
                if ~all(size(inter.coupling.eigs)==size(inter.coupling.euler))
                    error('validate: inter.coupling.eigs and inter.coupling.euler arrays must have the same dimension.');
                end
                
                % Make sure all non-empty elements are 3-vectors
                for n=1:spin_system.comp.nspins
                    for k=1:spin_system.comp.nspins
                
                    % For eigenvalues
                    if ((~all(size(inter.coupling.eigs{n,k})==[1 3]))&&(~all(size(inter.coupling.eigs{n,k})==[0 0])))||(~isnumeric(inter.coupling.eigs{n,k}))
                        error('validate: non-empty elements of inter.coupling.eigs must be 1x3 vectors.');
                    end
                    
                    % For Euler angles
                    if ((~all(size(inter.coupling.euler{n,k})==[1 3]))&&(~all(size(inter.coupling.euler{n,k})==[0 0])))||(~isnumeric(inter.coupling.euler{n,k}))
                        error('validate: non-empty elements of inter.coupling.euler must be 1x3 vectors.');
                    end
                    
                    end
                end
                
            end
            
            % Check the matrix specification
            if isfield(inter.coupling,'matrix')
                
                % Check the type
                if ~iscell(inter.coupling.matrix)
                    error('validate: inter.coupling.matrix must be a cell array of 3x3 matrices.');
                end
                
                % Check the dimensions
                if ~all(size(inter.coupling.matrix)==[spin_system.comp.nspins spin_system.comp.nspins])
                    error('validate: both dimensions of inter.coupling.matrix array should match the number of spins.');
                end
                
                % Make sure all non-empty elements are 3x3 matrices
                for n=1:spin_system.comp.nspins
                    for k=1:spin_system.comp.nspins
                        if ((~all(size(inter.coupling.matrix{n,k})==[3 3]))&&(~all(size(inter.coupling.matrix{n,k})==[0 0])))||(~isnumeric(inter.coupling.matrix{n,k}))
                            error('validate: non-empty elements of inter.coupling.matrix must be 3x3 matrices.');
                        end
                    end
                end
                
            end
            
            % Check the scalars
            if isfield(inter.coupling,'scalar')
                
                 % Check the type
                if ~iscell(inter.coupling.scalar)
                    error('validate: inter.coupling.scalar must be a cell array of numbers.');
                end
                
                % Check the dimensions
                if ~all(size(inter.coupling.scalar)==[spin_system.comp.nspins spin_system.comp.nspins])
                    error('validate: both dimensions of inter.coupling.scalar array should match the number of spins.');
                end
                
                % Make sure all non-empty elements are numbers
                for n=1:spin_system.comp.nspins
                    for k=1:spin_system.comp.nspins
                        if ((numel(inter.coupling.scalar{n,k})~=1)&&(numel(inter.coupling.scalar{n,k})~=0))||(~isnumeric(inter.coupling.scalar{n,k}))
                            error('validate: non-empty elements of inter.coupling.scalar must be numbers.');
                        end
                    end
                end
                
            end
                        
        end
        
    case 'create_relaxation_input'           % Done
        
        % Move to the natural name
        inter=parameter;
        
        % Check temperature - negative and complex values allowed
        if isfield(inter,'temperature')
            
            % Check type and dimension
            if (~isnumeric(inter.temperature))||(numel(inter.temperature)~=1)
                error('validate: inter.tempearture must be a number.');
            end
            
        end
        
        % Check relaxation theory
        if isfield(inter,'relaxation')
            
            % Check type
            if ~ischar(inter.relaxation)
                error('validate: inter.relaxation must be a string.');
            end
            
            % Check contents
            if ~any(strcmp(inter.relaxation,{'none','damp','t1_t2','redfield','sle'}))
                error('validate: allowed values for inter.relaxation are `none`, `damp`, `t1_t2`, `redfield` and `sle`.');
            end
                        
            % Enforce no relaxation in solids
            if (~strcmp(spin_system.inter.regime,'liquid'))&&strcmp(inter.relaxation,'redfield')
                error('validate: Redfield theory is not applicable to solid state.');
            end
            
            % Enforce correlation time with Redfield theory
            if strcmp(inter.relaxation,'redfield')&&(~isfield(inter,'tau_c'))
                error('validate: correlation time must be specified with Redfield theory.');
            end
            
            % Enforce correlation time with SLE
            if strcmp(inter.relaxation,'sle')&&(~isfield(inter,'tau_c'))
                error('validate: correlation time must be specified with SLE.');
            end
            
            % Enforce damping rate with non-selective damping
            if strcmp(inter.relaxation,'damp')&&(~isfield(inter,'damp_rate'))
                error('validate: damping rate must be specified with non-selective damping.');
            end
            
            % Enforce R1 and R2 rates with extended T1,T2 theory
            if strcmp(inter.relaxation,'t1_t2')&&((~isfield(inter,'r1_rates'))||(~isfield(inter,'r2_rates')))
                error('validate: R1 and R2 rates must be specified with extended T1,T2 relaxation theory.');
            end
            
        end
        
        % Check correlation time
        if isfield(inter,'tau_c')
            
            % Check type and dimension
            if (~isnumeric(inter.tau_c))||(numel(inter.tau_c)>3)||(numel(inter.tau_c)==0)
                error('validate: inter.tau_c must be a vector of size 1, 2 or 3.');
            end
            
            % Check value
            if (~isreal(inter.tau_c))||any(inter.tau_c<0)
                error('validate: inter.tau_c must have non-negative real elements.');
            end
            
            % Enforce Redfield theory if tau_c is specified
            if (~isfield(inter,'relaxation'))||(~any(strcmp(inter.relaxation,{'redfield','sle'})))
                error('validate: inter.tau_c requires Redfield theory or SLE.');
            end
            
        end
        
        % Check term retention
        if isfield(inter,'rlx_keep')
            
            % Check type
            if ~ischar(inter.rlx_keep)
                error('validate: inter.rlx_keep must be a string.');
            end
            
            % Check contents
            if ~any(strcmp(inter.rlx_keep,{'diagonal','kite','secular','full'}))
                error('validate: allowed values for inter.rlx_keep are `diagonal`, `kite`, `secular` and `full`.');
            end
                        
        end
        
        % Check DFS retention
        if isfield(inter,'rlx_dfs')
            
            % Check type
            if ~ischar(inter.rlx_dfs)
                error('validate: inter.rlx_dfs must be a string.');
            end
            
            % Check contents
            if ~any(strcmp(inter.rlx_dfs,{'keep','ignore'}))
                error('validate: allowed values for inter.rlx_dfs are `keep` and `ignore`.');
            end
                        
        end
        
        % Check equilibrium switch
        if isfield(inter,'equilibrium')
            
            % Check type
            if ~ischar(inter.equilibrium)
                error('validate: inter.equilibrium must be a string.');
            end
            
            % Check contents
            if ~any(strcmp(inter.equilibrium,{'zero','thermal'}))
                error('validate: allowed values for inter.equilibrium are `zero` and `thermal`.');
            end
                        
        end
        
        % Check non-selective damping rate
        if isfield(inter,'damp_rate')
            
            % Check type and dimension
            if (~isnumeric(inter.damp_rate))||(numel(inter.damp_rate)~=1)
                error('validate: inter.damp_rate must be a number.');
            end
            
            % Check value
            if (~isreal(inter.damp_rate))||(inter.damp_rate<0)
                error('validate: inter.damp_rate must be a non-negative real number.');
            end
            
            % Enforce non-selective damping if tau_c is specified
            if ~strcmp(inter.relaxation,'damp')
                error('validate: damp rate can only be specified with non-selective damping.');
            end
            
        end
        
        % Check the R1 rates
        if isfield(inter,'r1_rates')
            
            % Check type and values
            if (~isnumeric(inter.r1_rates))||any(inter.r1_rates<0)||any(~isreal(inter.r1_rates))
                error('validate: inter.r1_rates must be an array of positive real numbers.');
            end
            
            % Check dimension
            if numel(inter.r1_rates)~=spin_system.comp.nspins
                error('validate: the number of elements in inter.r1_rates must be equal to the number of spins.');
            end
            
            % Enforce T1,T2 theory if R1 rates are specified
            if ~strcmp(inter.relaxation,'t1_t2')
                error('validate: R1 rates can only be specified with extended T1,T2 relaxation theory.');
            end
            
        end
        
        % Check the R2 rates
        if isfield(inter,'r2_rates')
            
            % Check type and values
            if (~isnumeric(inter.r2_rates))||any(inter.r2_rates<0)||any(~isreal(inter.r2_rates))
                error('validate: inter.r2_rates must be an array of positive real numbers.');
            end
            
            % Check dimension
            if numel(inter.r2_rates)~=spin_system.comp.nspins
                error('validate: the number of elements in inter.r2_rates must be equal to the number of spins.');
            end
            
            % Enforce T1,T2 theory if R2 rates are specified
            if ~strcmp(inter.relaxation,'t1_t2')
                error('validate: R2 rates can only be specified with extended T1,T2 relaxation theory.');
            end
            
        end
        
    case 'create_chemistry_input'            % Done
        
        % Move to the natural name
        inter=parameter;
        
        % Check the chemical reactions
        if isfield(inter,'chem')
            
            % Check exchange matrix
            if isfield(inter.chem,'exchange')
                if ~all(size(inter.chem.exchange)==[spin_system.comp.nspins spin_system.comp.nspins])
                    error('validate: both dimensions of inter.chem.exchange matrix must be equal to the number of spins.');
                end
            end
            
        end
        
    case 'create_symmetry_input'             % Done
        
        % Move to the natural name
        sys=parameter;
        
        % Check symmetry parameters
        if isfield(sys,'sym_group')
            
            % Check the type
            if ~iscell(sys.sym_group)
                error('validate: sys.sym_group must be a cell array of strings.');
            end
            
            % Check that sys.sym_spins exists
            if ~isfield(sys,'sym_spins')
                error('validate: sys.sym_spins must be specified alongside sys.sym_group.');
            end
            
            % Check the type
            if ~iscell(sys.sym_spins)
                error('validate: sys.sym_spins must be a cell array of vectors.');
            end
            
            % Check the dimensions
            if numel(sys.sym_spins)~=numel(sys.sym_group)
                error('validate: sys.sym_group and sys.sym_spins arrays must have the same number of elements.');
            end
            
            % Check the spin indices
            for m=1:length(sys.sym_spins)
                
                % Check for sense
                if any(sys.sym_spins{m}>spin_system.comp.nspins)||any(sys.sym_spins{m}<1)
                    error('validate: incorrect spin labels in sys.sym_spins.');
                end
                
                % Check for intersections
                for n=1:length(sys.sym_spins)
                    if (n~=m)&&(~isempty(intersect(sys.sym_spins{m},sys.sym_spins{n})))
                        error('validate: same spin is listed in multiple symmetry groups in sys.sym_spins.');
                    end
                end
                
            end
            
            % Check the group names
            for n=1:length(sys.sym_group)
                if ~any(strcmp(sys.sym_group{n},{'S2','Ci','C2','Cs','S3','C3v','D3','S4','Td','S5','S6','D2h'}))
                    error('validate: the group requested in sys.sym_group is not available.');
                end
            end
            
            % Check the irrep switch
            if isfield(sys,'sym_a1g_only')&&(~isnumeric(sys.sym_a1g_only))&&(~islogical(sys.sym_a1g_only))&&((sys.sym_a1g_only~=1)||(sys.sym_a1g_only~=0))
                error('validate: the allowed values for sys.sym_a1g_only are 0 and 1.');
            end
            
        else
            
            % Enforce no sys.sym_spins without sys.sym_group
            if isfield(sys,'sym_spins')
                error('validate: sys.sym_group must be specified alongside sys.sym_spins.');
            end
            
            % Enforce no sys.sym_a1g_only without sys.sym_group
            if isfield(sys,'sym_a1g_only')
                error('validate: sys.sym_group must be specified alongside sys.sym_a1g_only.');
            end
            
        end
        
    case 'basis_input'                       % Done
        
        % Move to the natural name
        bas=parameter;
        
        % Check bas.mode
        if ~isfield(bas,'mode')
            error('validate: basis specification in bas.mode is required.');
        elseif ~ischar(bas.mode)
            error('validate: bas.mode must be a string.');
        elseif ~any(strcmp(bas.mode,{'complete','IK-0','IK-1','IK-2','ESR-1','ESR-2'}))
            error('validate: the basis set requested in bas.mode is not available.');
        end
        
        % Check bas.level
        if isfield(bas,'level')
            if (~isnumeric(bas.level))||(round(bas.level)~=bas.level)||(bas.level<1)
                error('validate: bas.level must be a positive integer.');
            end
        end
        
        % Check bas.space_level
        if isfield(bas,'space_level')
            if  (~isnumeric(bas.space_level))||(round(bas.space_level)~=bas.space_level)||(bas.space_level<1)
                error('validate: bas.space_level must be a positive integer.');
            end
        end
        
        % Check bas.manual
        if isfield(bas,'manual')
            if ~(islogical(bas.manual)||isnumeric(bas.manual))
                error('validate: bas.manual must be a logical matrix.');
            elseif size(bas.manual,2)~=spin_system.comp.nspins
                error('validate: the number of rows in bas.manual must be equal to the number of spins in the system.')
            end
        end
        
        % Check bas.cnes
        if isfield(bas,'cnes')
            if (~(isnumeric(bas.cnes)||islogical(bas.cnes)))||(bas.cnes~=0&&bas.cnes~=1)
                error('validate: bas.cnes must be 1 or 0.');
            end
        end
        
    case 'isotropic_couplings'
        
        % Move to the natural name
        coupling_iso=parameter;
        
        % Make sure all couplings are real
        for n=1:size(coupling_iso,1)
            for k=1:size(coupling_iso,2)
                if (numel(coupling_iso{n,k})>0)&&(abs(imag(coupling_iso{n,k}))>max_imag)
                    error('KERNEL INTEGRITY ERROR: significant imaginary component detected in a scalar coupling.');
                end
            end
        end
        
    case 'hermitian_superoperator'
        
        % Move to the natural name
        H=parameter;
        
        % Check that the matrix is hermitian
        if max(max(abs(H-H')))>max_nonherm
            error('KERNEL INTEGRITY ERROR: a superoperator that is meant to be Hermitian is not.');
        end
        
    otherwise
        
        % Do nothing
        
end

% Potentially, a government is the most dangerous threat to man's rights:
% it holds a legal monopoly on the use of physical force against legally
% disarmed victims.
%
% Ayn Rand















