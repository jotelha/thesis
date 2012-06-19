% Parameter and option preprocessing. See the manual for the input options.
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk
% kpervushin@ntu.edu.sg
% luke.edwards@chem.ox.ac.uk

function spin_system=create(sys,inter)

% Bomb out if Matlab is too old
if verLessThan('matlab', '7.12')
    error('create: your Matlab is too old to run Spinach - version 2011a or later required.');
end

%% Runtime settings and tolerances

% Input validation
validate([],sys,'create_runtime_input');

% Output switches
if isfield(sys,'output')
    % Use ASCII log file if asked
    spin_system.sys.output=sys.output;
    if isfield(sys,'logfile')
        spin_system.sys.logfile=sys.logfile;
    end
else
    % Use Matlab console by default
    spin_system.sys.output='console';
end

% Paths and locations
root_dir=which('h_superop');
if isempty(root_dir)
    error('create: Spinach directory not found - please follow installation instructions.');
else
    spin_system.sys.root_dir=root_dir(1:end-19);
    report(spin_system,['create: Spinach root directory determined to be ' spin_system.sys.root_dir]);
end

% Internal algorithms to disable
if isfield(sys,'disable')
    spin_system.sys.disable=sys.disable;
else
    spin_system.sys.disable={};
end

% Architecture-specific slashes
if ispc
    spin_system.sys.slash='\';
elseif isunix
    spin_system.sys.slash='/';
end

% Spinach version banner
banner(spin_system,'version_banner');

% Cut-offs and tolerances
spin_system=tolerances(spin_system,sys);

% Spin system banner
banner(spin_system,'spin_system_banner');

%% Isotopes, multiplicities and text labels

% Validate the input
validate(spin_system,sys,'create_comp_input');

% Absorb the operator caching identifier
if isfield(sys,'identifier')
    spin_system.sys.identifier=sys.identifier;
    report(spin_system,['create: operator caching enabled using "' spin_system.sys.identifier '" identifier.']);
end

% Absorb the liquid/crystal/powder regime switch
if isfield(sys,'regime')
    spin_system.inter.regime=sys.regime;
else
    spin_system.inter.regime='liquid';
end

% Inform the user about the connectivity analysis to be run
switch spin_system.inter.regime
    case 'liquid'
        report(spin_system,'create: LIQUID REGIME, spin system connectivity will be inferred from scalar couplings.');
    case {'crystal','powder'}
        report(spin_system,'create: SOLID REGIME, full coupling tensors will be used to infer spin system connectivity.');
    otherwise
        error('create: unknown regime.');
end

% Number and types of spins
spin_system.comp.isotopes=sys.isotopes;
spin_system.comp.nspins=numel(spin_system.comp.isotopes);
report(spin_system,['create: ' num2str(spin_system.comp.nspins) ' spins in the simulation.']);

% Text labels for spins
if isfield(sys,'labels')
    spin_system.comp.labels=sys.labels;
else
    spin_system.comp.labels=cell(spin_system.comp.nspins,1);
end

% Multiplicities and magnetogyric ratios
spin_system.comp.mults=zeros(1,spin_system.comp.nspins);
spin_system.comp.gamma=zeros(1,spin_system.comp.nspins);
for n=1:spin_system.comp.nspins
    [spin_system.comp.gamma(n),spin_system.comp.mults(n)]=spin(sys.isotopes{n});
end

% Order matrix
if isfield(inter,'order_matrix')
    spin_system.inter.order_matrix=inter.order_matrix;
else
    spin_system.inter.order_matrix=zeros(3);
end

%% Zeeman interactions

% Validate the input
validate(spin_system,{sys,inter},'create_zeeman_input');

% B0 magnetic induction
if isfield(sys,'magnet')
    spin_system.inter.magnet=sys.magnet;
    report(spin_system,['create: magnetic induction of ' num2str(spin_system.inter.magnet,'%0.5g') ' Tesla ('...
                        num2str(1e-6*spin_system.inter.magnet*spin('1H')/(2*pi),'%0.5g') ' MHz proton frequency, '...
                        num2str(1e-9*spin_system.inter.magnet*spin('E')/(2*pi),'%0.5g') ' GHz electron frequency).']);
else
    report(spin_system,'create: WARNING - magnet induction not specified - zero assumed.');
    spin_system.inter.magnet=0;
end

% Base frequencies
spin_system.inter.basefrq=spin_system.comp.gamma*spin_system.inter.magnet;

% Preallocate Zeeman tensor array
spin_system.inter.zeeman.matrix=mat2cell(zeros(3*spin_system.comp.nspins,3),3*ones(spin_system.comp.nspins,1));

% Process Zeeman interactions
if isfield(inter,'zeeman')

    % Absorb eigenvalues and Euler angles
    if isfield(inter.zeeman,'eigs')
        for n=1:spin_system.comp.nspins
            if significant(inter.zeeman.eigs{n},'vector',0)
                if (~isfield(inter.zeeman,'euler'))||isempty(inter.zeeman.euler{n})
                    S=eye(3);
                else
                    S=euler2dcm(inter.zeeman.euler{n});
                end
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+S*diag(inter.zeeman.eigs{n})*S';
            end
        end
    end
    
    % Absorb tensors
    if isfield(inter.zeeman,'matrix')
        for n=1:spin_system.comp.nspins
            if significant(inter.zeeman.matrix{n},'tensor',0)
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+inter.zeeman.matrix{n};
            end
        end
    end
    
    % Absorb scalars
    if isfield(inter.zeeman,'scalar')
         for n=1:spin_system.comp.nspins
            if significant(inter.zeeman.scalar{n},'scalar',0)
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+eye(3)*inter.zeeman.scalar{n};
            end
        end
    end

    % Report back to the user
    summary(spin_system,'zeeman','create: summary of non-zero Zeeman interactions (ppm for nuclei, g-tensor for electrons)');

    % Convert to angular frequencies
    for n=1:spin_system.comp.nspins
        switch spin_system.comp.isotopes{n}(1)
            case 'E'
                % For electrons, assume that the g-factor is given and compute the offset from the free electron g-factor
                spin_system.inter.zeeman.matrix{n}=(spin_system.inter.zeeman.matrix{n}-eye(3)*spin_system.tols.freeg)*(spin_system.inter.basefrq(n)/spin_system.tols.freeg);
            otherwise
                % For nuclei, assume that the chemical shift is given and compute the corresponding offset 
                spin_system.inter.zeeman.matrix{n}=(1e-6)*spin_system.inter.zeeman.matrix{n}*spin_system.inter.basefrq(n);
        end
    end
    
    % Symmetrize the result
    for n=1:spin_system.comp.nspins
        spin_system.inter.zeeman.matrix{n}=symmetrize(spin_system,spin_system.inter.zeeman.matrix{n});
    end

    % Clean up the result
    for n=1:spin_system.comp.nspins
        if negligible(spin_system.inter.zeeman.matrix{n},'tensor',spin_system.tols.inter_cutoff)
            spin_system.inter.zeeman.matrix{n}=zeros(3);
        end
    end

    % Report back to the user
    summary(spin_system,'zeeman','create: summary of non-zero Zeeman interactions (angular frequencies)');
    
else
    
    % Warn the user that Zeeman tensors have not been found
    report(spin_system,'create: WARNING - no Zeeman interactions supplied, magnet frequencies assumed.');
    
end

%% Spin-spin couplings

% Validate the input
validate(spin_system,inter,'create_coupling_input');

% Preallocate the tensor array
spin_system.inter.coupling.matrix=mat2cell(zeros(3*spin_system.comp.nspins),3*ones(spin_system.comp.nspins,1),3*ones(spin_system.comp.nspins,1));

% Process the coordinates
if isfield(inter,'coordinates')

    % Absorb the coordinates
    spin_system.inter.coordinates=inter.coordinates;
    
    % Preallocate the arrays
    spin_system.inter.distvectors=cell(spin_system.comp.nspins);
    spin_system.inter.distmatrix=cell(spin_system.comp.nspins);
    spin_system.inter.proxmatrix=zeros(spin_system.comp.nspins);
    
    % Loop over pairs of spins
    for n=1:spin_system.comp.nspins
        for k=1:spin_system.comp.nspins
            if (~isempty(spin_system.inter.coordinates{n}))&&(~isempty(spin_system.inter.coordinates{k}))&&(n~=k)
                
                % Write the distance vectors
                spin_system.inter.distvectors{n,k}=spin_system.inter.coordinates{n}-spin_system.inter.coordinates{k};
                
                % Write the distances
                spin_system.inter.distmatrix{n,k}=norm(spin_system.inter.distvectors{n,k});
                
                % Write the coupling matrices
                if spin_system.inter.distmatrix{n,k}<spin_system.tols.prox_cutoff
                    
                    % Write the proximity matrix
                    spin_system.inter.proxmatrix(n,k)=1;
                    
                    % Get the distance ort
                    ort=spin_system.inter.distvectors{n,k}/norm(spin_system.inter.distvectors{n,k});
                    
                    % Compute the dipolar interaction prefactor (0.5 due to double counting)
                    A=0.5*spin_system.comp.gamma(n)*spin_system.comp.gamma(k)*1.054571628e-34*1e-7/(norm(spin_system.inter.distvectors{n,k})*1e-10)^3;
                    
                    % Compute the dipolar coupling matrix
                    D=A*[1-3*ort(1)*ort(1)   -3*ort(1)*ort(2)   -3*ort(1)*ort(3);
                          -3*ort(2)*ort(1)  1-3*ort(2)*ort(2)   -3*ort(2)*ort(3);
                          -3*ort(3)*ort(1)   -3*ort(3)*ort(2)  1-3*ort(3)*ort(3)];
                      
                    % Add to the total
                    spin_system.inter.coupling.matrix{n,k}=spin_system.inter.coupling.matrix{n,k}+symmetrize(spin_system,D);
                    
                end
                
            end
        end
    end
    
    % Report back to the user
    summary(spin_system,'coordinates','create: coordinates (Angstrom)');
    summary(spin_system,'distances','create: distance matrix (Angstrom)');
    summary(spin_system,'dipolar','create:           summary of dipolar interactions     ');

else
    
    % Warn the user that the coordinates have not been found
    report(spin_system,'create: WARNING - no coordinates given, point dipolar interactions assumed to be zero.');
    spin_system.inter.proxmatrix=zeros(spin_system.comp.nspins);
    
end

% Absorb the user-specified couplings
if isfield(inter,'coupling')
    
    % Absorb eigenvalues and Euler angles
    if isfield(inter.coupling,'eigs')
        for n=1:spin_system.comp.nspins
            for k=1:spin_system.comp.nspins
                if significant(inter.coupling.eigs{n,k},'vector',0)
                    if (~isfield(inter.coupling,'euler'))||isempty(inter.coupling.euler{n,k})
                        S=eye(3);
                    else
                        S=euler2dcm(inter.coupling.euler{n,k});
                    end
                    spin_system.inter.coupling.matrix{n,k}=spin_system.inter.coupling.matrix{n,k}+2*pi*S*diag(inter.coupling.eigs{n,k})*S';
                end
            end
        end
    end
 
    % Absorb coupling tensors
    if isfield(inter.coupling,'matrix')
        for n=1:spin_system.comp.nspins
            for k=1:spin_system.comp.nspins
                if significant(inter.coupling.matrix{n,k},'tensor',0)
                    spin_system.inter.coupling.matrix{n,k}=spin_system.inter.coupling.matrix{n,k}+2*pi*inter.coupling.matrix{n,k};
                end
            end
        end
    end
        
    % Absorb scalar couplings
    if isfield(inter.coupling,'scalar')
        for n=1:spin_system.comp.nspins
            for k=1:spin_system.comp.nspins
                if significant(inter.coupling.scalar{n,k},'scalar',0)
                    spin_system.inter.coupling.matrix{n,k}=spin_system.inter.coupling.matrix{n,k}+2*pi*eye(3)*inter.coupling.scalar{n,k};
                end
            end
        end
    end
    
else
    
    % Warn the user that the couplings have not been found
    report(spin_system,'create: WARNING - no couplings given, zeros assumed.');

end

% Order up the coupling tensors
for n=2:spin_system.comp.nspins
    for k=1:(n-1)
        if ~isempty(spin_system.inter.coupling.matrix{n,k})
            if ~isempty(spin_system.inter.coupling.matrix{k,n})
                spin_system.inter.coupling.matrix{k,n}=spin_system.inter.coupling.matrix{k,n}+spin_system.inter.coupling.matrix{n,k};
            else
                spin_system.inter.coupling.matrix{k,n}=spin_system.inter.coupling.matrix{n,k};
            end
            spin_system.inter.coupling.matrix{n,k}=zeros(3);
        end
    end
end

% Symmetrize the coupling tensors
for n=1:spin_system.comp.nspins
    for k=1:spin_system.comp.nspins
        spin_system.inter.coupling.matrix{n,k}=symmetrize(spin_system,spin_system.inter.coupling.matrix{n,k});
    end
end

% Clean up the result
for n=1:spin_system.comp.nspins
    for k=1:spin_system.comp.nspins
        if negligible(spin_system.inter.coupling.matrix{n,k},'tensor',spin_system.tols.inter_cutoff)
            spin_system.inter.coupling.matrix{n,k}=[];
        end
    end
end

% Report back to the user
summary(spin_system,'couplings','create: summary of spin-spin couplings (angular frequencies)');

%% Relaxation parameters

% Validate the input
validate(spin_system,inter,'create_relaxation_input');

% Temperature
if ~isfield(inter,'temperature')
    inter.temperature=0;
end
spin_system.rlx.temperature=inter.temperature;
report(spin_system,['create: spin temperature: ' num2str(spin_system.rlx.temperature) ' Kelvin.']);

% Relaxation superoperator
if isfield(inter,'relaxation')
    spin_system.rlx.theory=inter.relaxation;
else
    spin_system.rlx.theory='none';
end
report(spin_system,['create: relaxation theory: ' spin_system.rlx.theory '.']);

% Rotational correlation times and rotational diffusion constants
if isfield(inter,'tau_c')
    spin_system.rlx.tau_c=inter.tau_c;
else
    spin_system.rlx.tau_c=0;
end
report(spin_system,['create: rotational correlation time(s): ' num2str(spin_system.rlx.tau_c) ' seconds.']);

% Terms to keep in the relaxation superoperator
if isfield(inter,'rlx_keep')
    spin_system.rlx.keep=inter.rlx_keep;
else
    spin_system.rlx.keep='kite';
end
report(spin_system,['create: terms to keep in the relaxation superoperator: ' spin_system.rlx.keep '.']);

% The fate of the dynamic frequency shift
if isfield(inter,'rlx_dfs')
    spin_system.rlx.dfs=inter.rlx_dfs;
else
    spin_system.rlx.dfs='ignore';
end
report(spin_system,['create: action to take on dynamic frequency shifts: ' spin_system.rlx.dfs '.']);

% Equilibrium state
if isfield(inter,'equilibrium')
    spin_system.rlx.equilibrium=inter.equilibrium;
else
    spin_system.rlx.equilibrium='zero';
end
report(spin_system,['create: equilibrium state to relax towards: ' spin_system.rlx.equilibrium '.']);

% Unit state warning for relaxation to thermal equilibrium
if strcmp(spin_system.rlx.equilibrium,'thermal');
    report(spin_system,'create: WARNING - the system will relax to thermal equilibirum.');
    report(spin_system,'create: WARNING - identity state populations will be forced to zero.');
end

% Phenomenological damping
if isfield(inter,'damp_rate')
    spin_system.rlx.damp_rate=inter.damp_rate;
else
    spin_system.rlx.damp_rate=0;
end
report(spin_system,['create: non-specific relaxation rate: ' num2str(spin_system.rlx.damp_rate) ' Hz.']);

% User-supplied R1 relaxation rates
if isfield(inter,'r1_rates')
    spin_system.rlx.r1=inter.r1_rates;
else
    spin_system.rlx.r1=zeros(spin_system.comp.nspins,1);
end

% User-supplied R2 relaxation rates
if isfield(inter,'r2_rates')
    spin_system.rlx.r2=inter.r2_rates;
else
    spin_system.rlx.r2=zeros(spin_system.comp.nspins,1);
end

% Report relaxation rates back to the user
if isfield(inter,'r1_rates')||isfield(inter,'r2_rates')
    summary(spin_system,'rlx_rates','create: summary of empirical relaxation rates (Hz)');
end

%% Spin system connectivity

% Connectivity matrix
switch spin_system.inter.regime
    % In liquids, use scalar couplings
    case 'liquid'
        spin_system.inter.conmatrix=sparse(abs(cellfun(@trace,spin_system.inter.coupling.matrix)/3)>spin_system.tols.inter_cutoff);
    % In solids, use full coupling tensors
    case {'crystal','powder'}
        spin_system.inter.conmatrix=sparse(cellfun(@norm,spin_system.inter.coupling.matrix)>spin_system.tols.inter_cutoff);
end

% Make sure each spin is connected and proximate to itself
spin_system.inter.conmatrix=double(spin_system.inter.conmatrix|eye(size(spin_system.inter.conmatrix)));
spin_system.inter.proxmatrix=double(spin_system.inter.proxmatrix|eye(size(spin_system.inter.proxmatrix)));

% Make sure connectivity and proximity are reciprocal
spin_system.inter.conmatrix=double(logical(spin_system.inter.conmatrix+spin_system.inter.conmatrix'));
spin_system.inter.proxmatrix=double(logical(spin_system.inter.proxmatrix+spin_system.inter.proxmatrix'));

% Issue a report to the user
report(spin_system,['create: connectivity matrix density ' num2str(100*nnz(spin_system.inter.conmatrix)/numel(spin_system.inter.conmatrix)) '%']);
report(spin_system,['create: proximity matrix density ' num2str(100*nnz(spin_system.inter.proxmatrix)/numel(spin_system.inter.proxmatrix)) '%']);

% Determine the number of independent subsystems
n_subsystems=max(scomponents(spin_system.inter.conmatrix|spin_system.inter.proxmatrix));

% Print a notice to the user
if n_subsystems>1
    report(spin_system,['create: the system contains ' num2str(n_subsystems) ' non-interacting subsystems.']);
end

%% Chemical processes

% Validate the input
validate(spin_system,inter,'create_chemistry_input');

% Set the kinetics assumptions
if isfield(inter,'chem')&&isfield(inter.chem,'type')
    spin_system.chem.type=inter.chem.type;
else
    spin_system.chem.type='intermolecular';
end

% Absorb the reaction rates
if isfield(inter,'chem')&&isfield(inter.chem,'exchange')
    
    % Absorb the rates
    [spin_system.chem.to,...
     spin_system.chem.from,...
     spin_system.chem.rate]=find(inter.chem.exchange);
     spin_system.chem.matrix=inter.chem.exchange;
     
    % Report back to the user
    summary(spin_system,'chemistry','create: chemical exchage reaction rates');
   
else
    
    spin_system.chem.to=[];
    spin_system.chem.from=[];
    spin_system.chem.rate=[];
    spin_system.chem.matrix=[];
    
end

% Absorb the radical recombination parameters
if isfield(inter,'chem')&&isfield(inter.chem,'rp_theory')

    % Absorb the theory
    spin_system.chem.rp_theory=inter.chem.rp_theory;
    report(spin_system,['create: radical recombination theory: ' spin_system.chem.rp_theory]);

    % Absorb the spins
    spin_system.chem.rp_spins=inter.chem.rp_spins;
    report(spin_system,['create: recombining electrons at positions: ' num2str(spin_system.chem.rp_spins)]);

    % Absorb the rates
    spin_system.chem.rp_rates=inter.chem.rp_rates;
    report(spin_system,['create: singlet recombination rate: ' num2str(spin_system.chem.rp_rates(1)) ' Hz.']);
    report(spin_system,['create: triplet recombination rate: ' num2str(spin_system.chem.rp_rates(2)) ' Hz.']);

else
    
    spin_system.chem.rp_theory='off';
    spin_system.chem.rp_spins=[];
    spin_system.chem.rp_rates=[];
    
end

%% Symmetry

% Validate the input
validate(spin_system,sys,'create_symmetry_input');

% Symmetry group
if isfield(sys,'sym_group');
    spin_system.comp.sym.group=sys.sym_group;
else
    spin_system.comp.sym.group={};
end

% Symmetry-related spins
if isfield(sys,'sym_spins')
    spin_system.comp.sym.spins=sys.sym_spins;
else
    spin_system.comp.sym.spins={};
end

% Irreducible representation composition
if isfield(sys,'sym_a1g_only')
    spin_system.comp.sym.a1g_only=sys.sym_a1g_only;
else
    spin_system.comp.sym.a1g_only=true();
end

% Report back to the user
if ~isempty(spin_system.comp.sym.group)
    summary(spin_system,'symmetry','create: symmetry summary');
end

end

% "Those who beat their swords into plowshares will till the soil for those
% who have not."
%
% Benjamin Franklin

