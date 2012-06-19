% Symmetrizes interaction tensors.
%
% This is required to avoid numerical instabilities associated with small
% asymmetries during the diagonalization process.
%
% ilya.kuprov@oerc.ox.ac.uk

function tensor=symmetrize(spin_system,tensor)

% If the asymmetry is significant, make a fuss
if norm(tensor-tensor')>spin_system.tols.inter_sym
    spin_system=click(spin_system,'forward');
    report(spin_system,'symmetrize: WARNING - significant asymmetry detected in a coupling tensor.');
    report(spin_system,['symmetrize: WARNING - symmetric part norm: ' num2str(norm((tensor+tensor')/2))]);
    report(spin_system,['symmetrize: WARNING - anti-symmetric part norm: ' num2str(norm((tensor-tensor')/2))]);
    report(spin_system,'symmetrize: WARNING - the tensor has been symmetrized.');
end

% Symmetrize the tensor
tensor=(tensor+tensor')/2;

end

% "Smoking - NO HYDROGEN!"
%
% A safety warning on Anatole Abragam's door