% Commutation superoperators.
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function A=c_superop(spin_system,opspec)

% Validate the input
validate(spin_system,opspec,'operator_specification')

if nnz(opspec)==0

    % Return zero superoperator
    A=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,0);

else

    % Call product superoperators and take the difference
    A=p_superop(spin_system,opspec,'left')-p_superop(spin_system,opspec,'right');

end

% Validate the output
validate(spin_system,A,'superoperator')
     
end

% "Are you thinking that death and taxes are our only certainty, Mr.
% Rearden? Well, there's nothing I can do about the first, but if I lift
% the burden of the second, men might learn to see the connection between
% the two and what a longer, happier life they have the power to achieve.
% They might learn to hold, not death and taxes, but life and production
% as their two absolutes and as the base of their moral code."
%
% Ayn Rand, "Atlas Shrugged"

