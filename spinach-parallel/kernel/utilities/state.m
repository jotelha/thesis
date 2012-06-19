% Generates state vectors from their user-friendly descriptions. Arguments:
%	
%     states - a cell array of strings. Each string may be 'E' (identity 
%              state), 'Lz', 'L+', 'L-' or 'Tl,m' (higher spin states as
%              irreducible spherical tensors; l and m are both integers).
%
%      spins - 
%              EITHER a cell array of integers specifying the numbers of spins
%                     to be generated in states given in the 'states' argument
%                     (this produces the corresponding multi-spin state)
%              
%              OR     an isotope specification: '13C', '15N', 'all' (this
%                     produces a sum of single-spin states on all the spins
%                     of the specified isotope)
% Example:
%
%    LzSp=state(spin_system,{'Lz','L+'},{1,2});
%
% would return the LzS+ state with Lz on spin 1 and L+ on spin 2.
%
% Example:
%
%    Sum_Lz=state(spin_system,'Lz','13C');
%
% would return the sum of Lz states on all carbons in the system.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function rho=state(spin_system,states,spins)

% Generate the unit state
unit_state=sparse(1,1,1,spin_system.bas.nstates,1);

% Apply the left side product superoperator
rho=operator(spin_system,states,spins,'left')*unit_state;

end

% It worked.
%
% J. Robert Oppenheimer (after witnessing the first atomic detonation) 

