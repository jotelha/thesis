% Hamiltonian superoperator and its rotational decomposition for Liouville
% space simulations. Parameters:
%
%      operator_type - can be set to:
%                      
%                      'left'  - produces left side product superoperator
%                      'right' - produces right side product superoperator
%                      'comm'  - produces commutation superoperator (default)
%
% Outputs:
%                            H - the isotropic part of the superoperator
%                            Q - twenty-five superoperators giving the
%                                irreducible components of the anisotropic
%                                part of the superoperator.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function [H,Q]=h_superop(spin_system,operator_type)

%% Preliminaries

% Click forward for output
spin_system=click(spin_system,'forward');

% Set the default for the type
if nargin==1
    operator_type='comm';
end

% Check that secularity assumptions have been set
if ~isfield(spin_system.inter.zeeman,'strength')
    error(['h_superop: spin_system.inter.zeeman.strength has not\n'...
           'been set. Either run secularity() with the appropriate\n'...
           'setting or define spin_system.inter.zeeman.strength\n'...
           'manually.']);
elseif ~isfield(spin_system.inter.coupling,'strength')
    error(['h_superop: spin_system.inter.coupling.strength has not\n'...
           'been set. Either run secularity() with the appropriate\n'...
           'setting or define spin_system.inter.coupling.strength\n'...
           'manually.']);
end

% Preallocate the relevant arrays
H=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,10*spin_system.bas.nstates);
if nargout==2
    Q=cell(5); T=mat2cell(zeros(15,3),[3 3 3 3 3],3);
    for L=1:5
        for S=1:5
            Q{L,S}=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,10*spin_system.bas.nstates);
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
        Lz=operator(spin_system,'Lz',L,operator_type);

        % Write the isotropic part
        zeeman_iso=trace(spin_system.inter.zeeman.matrix{L})/3;
        if significant(zeeman_iso,'scalar',spin_system.tols.inter_cutoff)
            if strcmp(spin_system.inter.zeeman.strength{L},'full')
                report(spin_system,['h_superop: full isotropic Zeeman term for spin ' num2str(L) '...']);
            else
                report(spin_system,['h_superop: offset isotropic Zeeman term for spin ' num2str(L) '...']);
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
                    report(spin_system,['h_superop: full anisotropic Zeeman term for spin ' num2str(L) '...']);
                    Lm=operator(spin_system,'L-',L,operator_type);
                    T{1}=[];
                    T{2}=-0.5*Lm';
                    T{3}=sqrt(2/3)*Lz;
                    T{4}=+0.5*Lm;
                    T{5}=[];
                case 'secular'
                    report(spin_system,['h_superop: secular anisotropic Zeeman term for spin ' num2str(L) '...']);
                    T{1}=[]; T{2}=[];
                    T{3}=sqrt(2/3)*Lz;
                    T{4}=[]; T{5}=[];
                otherwise
                    error('h_superop: unknown Zeeman interaction strength specification.');
            end
            
            % Print the orientation functions
            report(spin_system,['           PHI( 2) ' num2str(phi_zeeman(1)/(2*pi),'%+0.5e') ' Hz']);
            report(spin_system,['           PHI( 1) ' num2str(phi_zeeman(2)/(2*pi),'%+0.5e') ' Hz']);
            report(spin_system,['           PHI( 0) ' num2str(phi_zeeman(3)/(2*pi),'%+0.5e') ' Hz']);
            report(spin_system,['           PHI(-1) ' num2str(phi_zeeman(4)/(2*pi),'%+0.5e') ' Hz']);
            report(spin_system,['           PHI(-2) ' num2str(phi_zeeman(5)/(2*pi),'%+0.5e') ' Hz']);
            
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
                    error('h_superop: quadratic couplings cannot have an isotropic part.'); 
                end
                
                LzSz=operator(spin_system,{'Lz','Lz'},{L,S},operator_type);
                switch spin_system.inter.coupling.strength{L,S}
                    case {'strong','secular'}
                        report(spin_system,['h_superop: strong isotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                        report(spin_system,['           (LxSx+LySy+LzSz) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                        LpSm=operator(spin_system,{'L+','L-'},{L,S},operator_type);
                        H=H+coupling_iso*(LzSz+0.5*(LpSm+LpSm'));
                    case {'weak','L-secular','S-secular'}
                        report(spin_system,['h_superop: weak isotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                        report(spin_system,['           LzSz x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                        H=H+coupling_iso*LzSz;
                    otherwise
                        error('h_superop: unknown coupling strength specification.');
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
                            report(spin_system,['h_superop: strong quadratic coupling term for spin ' num2str(L) '...']);
                            T{1}=operator(spin_system,'T2,+2',L,operator_type);
                            T{2}=operator(spin_system,'T2,+1',L,operator_type);
                            T{3}=operator(spin_system,'T2,0',L,operator_type);
                            T{4}=operator(spin_system,'T2,-1',L,operator_type);
                            T{5}=operator(spin_system,'T2,-2',L,operator_type);
                        case 'secular'
                            report(spin_system,['h_superop: secular quadratic coupling term for spin ' num2str(L) '...']);
                            T{1}=[]; T{2}=[];
                            T{3}=operator(spin_system,'T2,0',L,operator_type);
                            T{4}=[]; T{5}=[];
                        case 'weak'
                            error('h_superop: weak quadratic interactions do not exist.');
                        otherwise
                            error('h_superop: unknown coupling strength specification.');
                    end
                    
                end
                
                % Write the tensors for bilinear anisotropic couplings
                if (L~=S)&&significant(phi_coupling,'vector',spin_system.tols.inter_cutoff)
                    
                    % Write the irreducible spherical tensors
                    switch spin_system.inter.coupling.strength{L,S}
                        case 'strong'
                            report(spin_system,['h_superop: strong anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=operator(spin_system,{'Lz','Lz'},{L,S},operator_type);
                            LpSm=operator(spin_system,{'L+','L-'},{L,S},operator_type);
                            LmSm=operator(spin_system,{'L-','L-'},{L,S},operator_type);
                            LmSz=operator(spin_system,{'L-','Lz'},{L,S},operator_type);
                            LzSm=operator(spin_system,{'Lz','L-'},{L,S},operator_type);
                            T{1}=+0.5*LmSm';
                            T{2}=-0.5*(LzSm'+LmSz');
                            T{3}=sqrt(2/3)*(LzSz-0.25*(LpSm+LpSm'));
                            T{4}=+0.5*(LzSm+LmSz);
                            T{5}=+0.5*LmSm;
                        case 'L-secular'
                            report(spin_system,['h_superop: L-secular anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=operator(spin_system,{'Lz','Lz'},{L,S},operator_type);
                            LzSm=operator(spin_system,{'Lz','L-'},{L,S},operator_type);
                            T{1}=[];
                            T{2}=-0.5*LzSm';
                            T{3}=sqrt(2/3)*LzSz;
                            T{4}=+0.5*LzSm;
                            T{5}=[];
                        case 'S-secular'
                            report(spin_system,['h_superop: S-secular anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=operator(spin_system,{'Lz','Lz'},{L,S},operator_type);
                            LmSz=operator(spin_system,{'L-','Lz'},{L,S},operator_type);
                            T{1}=[];
                            T{2}=-0.5*LmSz';
                            T{3}=sqrt(2/3)*LzSz;
                            T{4}=+0.5*LmSz;
                            T{5}=[];
                        case 'secular'
                            report(spin_system,['h_superop: secular anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=operator(spin_system,{'Lz','Lz'},{L,S},operator_type);
                            LpSm=operator(spin_system,{'L+','L-'},{L,S},operator_type);
                            T{1}=[]; T{2}=[];
                            T{3}=sqrt(2/3)*(LzSz-0.25*(LpSm+LpSm'));
                            T{4}=[]; T{5}=[];
                        case 'weak'
                            report(spin_system,['h_superop: weak anisotropic coupling term for spins ' num2str(L) ',' num2str(S) '...']);
                            LzSz=operator(spin_system,{'Lz','Lz'},{L,S},operator_type);
                            T{1}=[]; T{2}=[];
                            T{3}=sqrt(2/3)*LzSz;
                            T{4}=[]; T{5}=[];
                        otherwise
                            error('h_superop: unknown coupling strength specification.');
                    end
                    
                end
                
                % Print the orientation functions
                if significant(phi_coupling,'vector',spin_system.tols.inter_cutoff)
                    report(spin_system,['           PHI( 2) ' num2str(phi_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           PHI( 1) ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           PHI( 0) ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           PHI(-1) ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           PHI(-2) ' num2str(phi_coupling(5)/(2*pi),'%+0.5e') ' Hz']);
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

% The soundtrack that accompanied much of the Spinach kernel programming is
% here: http://spindynamics.org/documents/soundtrack.rar
%
% The password is "superoperator".

