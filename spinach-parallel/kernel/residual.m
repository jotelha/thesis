% Sets up interaction tensors under partial ordering in a liquid crystal with
% the user-supplied order matrix.
%
% This module is only applicable to high-field NMR spectroscopy.
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@chem.ox.ac.uk

function spin_system=residual(spin_system)

% Click forward for output
spin_system=click(spin_system,'forward');

% Process Zeeman interactions
for n=1:spin_system.comp.nspins
    if significant(spin_system.inter.zeeman.matrix{n},'tensor',spin_system.tols.inter_cutoff)
        
        % Obtain isotropic part
        iso=eye(3)*trace(spin_system.inter.zeeman.matrix{n})/3;
        
        % Calculate residual order
        extra_zz=trace(spin_system.inter.order_matrix*(spin_system.inter.zeeman.matrix{n}-iso));
        
        % Update the Zeeman tensor
        spin_system.inter.zeeman.matrix{n}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3]);
        
    end
end

% Process the spin-spin couplings
for n=1:spin_system.comp.nspins
    for k=1:spin_system.comp.nspins
        if significant(spin_system.inter.coupling.matrix{n,k},'tensor',spin_system.tols.inter_cutoff)
            
            % Obtain isotropic part
            iso=trace(spin_system.inter.coupling.matrix{n,k})*eye(3)/3;
            
            % Calculate residual order
            extra_zz=trace(spin_system.inter.order_matrix*(spin_system.inter.coupling.matrix{n,k}-iso));
            
            % Update the coupling tensor
            spin_system.inter.coupling.matrix{n,k}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3]);
            
        end
    end
end

% Report back to the user
report(spin_system,'residual: all interaction anisotropies have been replaced by their residuals.');

% Click backward for output
spin_system=click(spin_system,'backward');

end

% If I'd observed all the rules, I'd never have got anywhere.
%
% Marilyn Monroe

