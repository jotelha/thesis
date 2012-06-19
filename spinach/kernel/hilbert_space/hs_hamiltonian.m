% Hamiltonian operator and its rotational decomposition for Hilbert space
% simulations.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function [H,Q]=hs_hamiltonian(spin_system)

%% Preliminaries

% Click forward for output
spin_system=click(spin_system,'forward');

% Check that secularity assumptions have been set
if ~isfield(spin_system.inter.zeeman,'strength')
    error('hs_hamiltonian:secularity',['hs_hamiltonian: spin_system.inter.zeeman.strength has not\n'...
          'been set. Either run secularity() with the appropriate\n'...
          'setting or define spin_system.inter.zeeman.strength\nmanually.'])
elseif ~isfield(spin_system.inter.coupling,'strength')
    error('hs_hamiltonian:secularity',['hs_hamiltonian: spin_system.inter.coupling.strength has not\n'...
          'been set. Either run secularity() with the appropriate\n'...
          'setting or define spin_system.inter.coupling.strength\nmanually.'])
end

% Preallocate the relevant arrays
dimension=prod(spin_system.comp.mults);
H=spalloc(dimension,dimension,10*dimension);
if nargout==2
    Q=cell(5); T=mat2cell(zeros(15,3),[3 3 3 3 3],3);
    for L=1:5
        for S=1:5
            Q{L,S}=spalloc(dimension,dimension,10*dimension);
        end
    end
end

%% Zeeman tensors

for L=1:spin_system.comp.nspins
    
    % Add carrier frequencies if necessary
    if strcmp(spin_system.inter.zeeman.strength{L},'full')
        spin_system.inter.zeeman.matrix{L}=spin_system.inter.zeeman.matrix{L}+eye(3)*spin_system.inter.basefrq(L);
    end
    
    if significant(spin_system.inter.zeeman.matrix{L},'tensor',spin_system.tols.inter_cutoff)
        
        % Compute the Lz operator
        Lz=hs_operator(spin_system,'Lz',L);
        
        % Write the isotropic part
        zeeman_iso=trace(spin_system.inter.zeeman.matrix{L})/3;
        if significant(zeeman_iso,'scalar',spin_system.tols.inter_cutoff)
            if strcmp(spin_system.inter.zeeman.strength{L},'full')
                report(spin_system,['hs_hamiltonian: full isotropic Zeeman term for spin ' num2str(L) '...']);
            else
                report(spin_system,['hs_hamiltonian: offset isotropic Zeeman term for spin ' num2str(L) '...']);
            end
            report(spin_system,['           Lz x ' num2str(zeeman_iso/(2*pi)) ' Hz']);
            H=H+zeeman_iso*Lz;
        end

        % Get the second rank spherical tensor components
        [~,~,phi_zeeman]=mat2sphten(spin_system.inter.zeeman.matrix{L});

        % Write the anisotropic part
        if (nargout==2)&&significant(phi_zeeman,'vector',spin_system.tols.inter_cutoff)
            
            % Write the irreducible spherical tensors
            switch spin_system.inter.zeeman.strength{L}
                case 'full'
                    report(spin_system,['hs_hamiltonian: full anisotropic Zeeman term for spin ' num2str(L) '...']);
                    Lm=hs_operator(spin_system,'L-',L);
                    T{1}=[];
                    T{2}=-0.5*Lm';
                    T{3}=sqrt(2/3)*Lz;
                    T{4}=+0.5*Lm;
                    T{5}=[];
                case 'secular'
                    report(spin_system,['hs_hamiltonian: secular anisotropic Zeeman term for spin ' num2str(L) '...']);
                    T{1}=[]; T{2}=[];
                    T{3}=sqrt(2/3)*Lz;
                    T{4}=[]; T{5}=[];
                otherwise
                    error('hs_hamiltonian: unknown Zeeman interaction strength specification.');
            end
            
            % Print the orientation functions
            report(spin_system,['           PHI(-2) ' num2str(phi_zeeman(1)/(2*pi)) ' Hz']);
            report(spin_system,['           PHI(-1) ' num2str(phi_zeeman(2)/(2*pi)) ' Hz']);
            report(spin_system,['           PHI( 0) ' num2str(phi_zeeman(3)/(2*pi)) ' Hz']);
            report(spin_system,['           PHI( 1) ' num2str(phi_zeeman(4)/(2*pi)) ' Hz']);
            report(spin_system,['           PHI( 2) ' num2str(phi_zeeman(5)/(2*pi)) ' Hz']);
            
            % Update the rotation basis
            for k=1:5
                if significant(T{k},'operator',eps)
                    for m=1:5
                        Q{k,m}=Q{k,m}+phi_zeeman(m)*T{k};
                    end
                end
            end
            
        end
        
    end
    
end

%% Coupling tensors

for L=1:spin_system.comp.nspins
    for S=1:spin_system.comp.nspins
        if significant(spin_system.inter.coupling.matrix{L,S},'tensor',spin_system.tols.inter_cutoff)
                       
            % Write the isotropic couplings
            coupling_iso=trace(spin_system.inter.coupling.matrix{L,S})/3;
            if significant(coupling_iso,'scalar',spin_system.tols.inter_cutoff)
                
                if L==S
                    error('hs_hamiltonian: quadratic couplings cannot have an isotropic part.'); 
                end
                
			    LzSz=hs_operator(spin_system,{'Lz','Lz'},{L,S});
                switch spin_system.inter.coupling.strength{L,S}
                    case {'strong','secular'}
                        report(spin_system,['hs_hamiltonian: strong isotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                        report(spin_system,['           (LxSx+LySy+LzSz) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                        LpSm=hs_operator(spin_system,{'L+','L-'},{L,S});
                        H=H+coupling_iso*(LzSz+0.5*(LpSm+LpSm'));
                    case {'weak','L-secular','S-secular'}
                        report(spin_system,['hs_hamiltonian: weak isotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                        report(spin_system,['           LzSz x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                        H=H+coupling_iso*LzSz;
                    otherwise
                        error('hs_hamiltonian: unknown coupling strength specification.');
                end
            end
            
            if nargout==2
                
                % Get the second rank spherical tensor components
                [~,~,phi_coupling]=mat2sphten(spin_system.inter.coupling.matrix{L,S});
                
                % Write the tensors for quadratic anisotropic couplings
                if (L==S)&&significant(phi_coupling,'vector',spin_system.tols.inter_cutoff)
                    
                    % Write the irreducible spherical tensors
                    switch spin_system.inter.coupling.strength{L,S}
                        case 'strong'
                            report(spin_system,['hs_hamiltonian: strong quadratic coupling term for spin ' num2str(L) '...']);
                            T{1}=hs_operator(spin_system,'T2,+2',L);
                            T{2}=hs_operator(spin_system,'T2,+1',L);
                            T{3}=hs_operator(spin_system,'T2,0',L);
                            T{4}=hs_operator(spin_system,'T2,-1',L);
                            T{5}=hs_operator(spin_system,'T2,-2',L);
                        case 'secular'
                            report(spin_system,['hs_hamiltonian: secular quadratic coupling term for spin ' num2str(L) '...']);
                            T{1}=[]; T{2}=[];
                            T{3}=sqrt(2/3)*hs_operator(spin_system,'T2,0',L);
                            T{4}=[]; T{5}=[];
                        case 'weak'
                            error('hs_hamiltonian: weak quadratic interactions do not exist.');
                        otherwise
                            error('hs_hamiltonian: unknown coupling strength specification.');
                    end
                    
                end
                
                % Write the tensors for bilinear anisotropic couplings
                if (L~=S)&&significant(phi_coupling,'vector',spin_system.tols.inter_cutoff)
                    
                    % Write the irreducible spherical tensors
                    switch spin_system.inter.coupling.strength{L,S}
                        case 'strong'
                            report(spin_system,['hs_hamiltonian: strong anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=hs_operator(spin_system,{'Lz','Lz'},{L,S});
                            LpSm=hs_operator(spin_system,{'L+','L-'},{L,S});
                            LmSm=hs_operator(spin_system,{'L-','L-'},{L,S});
                            LmSz=hs_operator(spin_system,{'L-','Lz'},{L,S});
                            LzSm=hs_operator(spin_system,{'Lz','L-'},{L,S});
                            T{1}=+0.5*LmSm';
                            T{2}=-0.5*(LzSm'+LmSz');
                            T{3}=sqrt(2/3)*(LzSz-0.25*(LpSm+LpSm'));
                            T{4}=+0.5*(LzSm+LmSz);
                            T{5}=+0.5*LmSm;
                        case 'L-secular'
                            report(spin_system,['hs_hamiltonian: L-secular anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=hs_operator(spin_system,{'Lz','Lz'},{L,S});
                            LzSm=hs_operator(spin_system,{'Lz','L-'},{L,S});
                            T{1}=[];
                            T{2}=-0.5*LzSm';
                            T{3}=sqrt(2/3)*LzSz;
                            T{4}=+0.5*LzSm;
                            T{5}=[];
                        case 'S-secular'
                            report(spin_system,['hs_hamiltonian: S-secular anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=hs_operator(spin_system,{'Lz','Lz'},{L,S});
                            LmSz=hs_operator(spin_system,{'L-','Lz'},{L,S});
                            T{1}=[];
                            T{2}=-0.5*LmSz';
                            T{3}=sqrt(2/3)*LzSz;
                            T{4}=+0.5*LmSz;
                            T{5}=[];
                        case 'secular'
                            report(spin_system,['h_superop: secular anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=hs_operator(spin_system,{'Lz','Lz'},{L,S});
                            LpSm=hs_operator(spin_system,{'L+','L-'},{L,S});
                            T{1}=[]; T{2}=[];
                            T{3}=sqrt(2/3)*(LzSz-0.25*(LpSm+LpSm'));
                            T{4}=[]; T{5}=[];
                        case 'weak'
                            report(spin_system,['h_superop: weak anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            T{1}=[]; T{2}=[];
                            T{3}=sqrt(2/3)*hs_operator(spin_system,{'Lz','Lz'},{L,S});
                            T{4}=[]; T{5}=[];
                        otherwise
                            error('h_superop: unknown coupling strength specification.');
                    end
                    
                end
                
                % Print the orientation functions
                if significant(phi_coupling,'vector',spin_system.tols.inter_cutoff)
                    report(spin_system,['           PHI(-2) ' num2str(phi_coupling(1)/(2*pi)) ' Hz']);
                    report(spin_system,['           PHI(-1) ' num2str(phi_coupling(2)/(2*pi)) ' Hz']);
                    report(spin_system,['           PHI( 0) ' num2str(phi_coupling(3)/(2*pi)) ' Hz']);
                    report(spin_system,['           PHI( 1) ' num2str(phi_coupling(4)/(2*pi)) ' Hz']);
                    report(spin_system,['           PHI( 2) ' num2str(phi_coupling(5)/(2*pi)) ' Hz']);
                end
                
                % Update the rotation basis
                if significant(phi_coupling,'vector',spin_system.tols.inter_cutoff)
                    for k=1:5
                        if significant(T{k},'operator',eps)
                            for m=1:5
                                Q{k,m}=Q{k,m}+phi_coupling(m)*T{k};
                            end
                        end
                    end
                end
                
            end
            
        end
    end
end

%% Miscellaneous things

% Kill the small terms
H=H.*(abs(H)>spin_system.tols.liouv_zero);
if nargout==2
    for n=1:5
        for k=1:5
            Q{n,k}=Q{n,k}.*(abs(Q{n,k})>spin_system.tols.liouv_zero);
        end
    end
end

end

% A cactus is a very disappointed cucumber.

