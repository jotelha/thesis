% Anticommutation superoperators.
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function A=a_superop(spin_system,opspec)

% Validate the input
validate(spin_system,opspec,'operator_specification')

% Call product superoperators and take the sum
A=p_superop(spin_system,opspec,'left')+p_superop(spin_system,opspec,'right');

% Validate the output
validate(spin_system,A,'superoperator')
     
end

% "There are 10 kinds of people in the world, those that understand
% ternary, those that don't, and those that confuse it with binary."

