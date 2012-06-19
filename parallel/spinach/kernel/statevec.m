% User-specified state vector. Converts a Spinach operator specification
% into a Liouville space state vector. Arguments:
%
%      opspec - Spinach operator specification (see the manual)
%      rho    - the resulting state vector
%
% If multiple lines are present in opspec, returns the sum of the corres-
% ponding state vectors.
%
% ilya.kuprov@oerc.ox.ac.uk

function rho=statevec(spin_system,opspec)

% Generate the unit state
unit_state=sparse(1,1,1,spin_system.bas.nstates,1);

% Apply the left side product superoperator
rho=p_superop(spin_system,opspec,'left')*unit_state;

end

% Evans boldly put 50 atm of ethylene in a cell with 25 atm of oxygen. The
% apparatus subsequently blew up, but luckily not before he had obtained
% the spectra shown in Figure 8.
%
% A.J. Mehrer and R.S. Mulliken, Chem. Rev. 69 (1969) 639-656