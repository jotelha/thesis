% Sets the interaction secularity assumptions for the coherent Liouvillian.
% Simulations of different experiments require different rotating frames
% and different secular approximations - those are stored in this function.
%
% ilya.kuprov@oerc.ox.ac.uk

function spin_system=secularity(spin_system,setting)

% Click forward for output
spin_system=click(spin_system,'forward');

% Report to the user
report(spin_system,['secularity: rotating frame assumptions set to "' setting '".']);

% Decide the secularity of Zeeman interactions
spin_system.inter.zeeman.strength=cell(1,spin_system.comp.nspins);
switch setting
    
    case {'nmr','esr'}
        
        % For high-field NMR and ESR, keep the offset ZZ term
        for n=1:spin_system.comp.nspins
            spin_system.inter.zeeman.strength{n}='secular';
        end
        
    case {'lowfield','keep_all'}
        
        % At low field or under direct instruction, keep everything
        for n=1:spin_system.comp.nspins
            spin_system.inter.zeeman.strength{n}='full';
        end
        
    case {'eseem','endor','solid_effect'}
        
        % For electrons, keep the offset ZZ term; for nuclei, keep everything
        for n=1:spin_system.comp.nspins
            if strcmp(spin_system.comp.isotopes{n}(1),'E')
                spin_system.inter.zeeman.strength{n}='secular';
            else
                spin_system.inter.zeeman.strength{n}='full';
            end
        end
        
    otherwise
        
        % With unknown settings, terminate
        error('secularity: unrecognised secularity setting.')
        
end

% Decide the secularity of couplings
spin_system.inter.coupling.strength=cell(spin_system.comp.nspins);
switch setting
    
    case {'nmr','esr'}
        
        % For high-field NMR and ESR, keep the secular terms for spins
        % belonging to the same isotope and weak terms otherwise.
        for n=1:spin_system.comp.nspins
            for k=1:spin_system.comp.nspins
                if strcmp(spin_system.comp.isotopes{n},spin_system.comp.isotopes{k})
                    spin_system.inter.coupling.strength{n,k}='secular';
                else
                    spin_system.inter.coupling.strength{n,k}='weak';
                end
            end
        end
        
    case {'lowfield','keep_all'}
        
        % At low field or under direct instruction, keep everything
        for n=1:spin_system.comp.nspins
            for k=1:spin_system.comp.nspins
                spin_system.inter.coupling.strength{n,k}='strong';
            end
        end
        
    case {'eseem','endor','solid_effect'}
        
        for n=1:spin_system.comp.nspins
            for k=1:spin_system.comp.nspins
                if strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    % Couplings from electron to nucleus should be L-secular
                    spin_system.inter.coupling.strength{n,k}='L-secular';
                elseif strcmp(spin_system.comp.isotopes{k}(1),'E')&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    % Couplings from nucleus to electron should be S-secular
                    spin_system.inter.coupling.strength{n,k}='S-secular';
                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    % Couplings between electrons should be secular
                    spin_system.inter.coupling.strength{n,k}='secular';
                else
                    % Couplings between nuclei should be strong
                    spin_system.inter.coupling.strength{n,k}='strong';
                end
            end
        end
        
    otherwise
        
        error('secularity: unrecognised `setting` argument.')
        
end

% Click backward for output
spin_system=click(spin_system,'backward');

end

% Humanity has advanced, when it has advanced, not because it has been sober,
% responsible, and cautious, but because it has been playful, rebellious, and immature.
%
% Tom Robbins

