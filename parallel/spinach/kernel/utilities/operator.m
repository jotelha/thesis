% Generates superoperators from their user-friendly descriptions. Arguments:
%	
%     operators - a cell array of strings. Each string may be 'E' (identity 
%                 operator), 'Lz', 'L+', 'L-' or 'Tl,m' (higher spin orders as
%                 irreducible spherical tensors; l and m are both integers).
%
%         spins - 
%                 EITHER a cell array of integers specifying the numbers of spins
%                        on which the operators given in the 'operators' argument
%                        operate (this produces the corresponding multi-spin
%                        superoperator)
%              
%                 OR     an isotope specification: '13C', '15N', 'all' (this
%                        produces a sum of single-spin operators on all the spins
%                        of the specified isotope)
%
% operator_type - can be set to:
%                      
%                 'left'  - produces left side product superoperator
%                 'right' - produces right side product superoperator
%                 'comm'  - produces commutation superoperator (default)
%
% Example:
%
%    LzSp=operator(spin_system,{'Lz','L+'},{1,2});
%
% would return the [LzS+, ] commutation superoperator with Lz on spin 1 and 
% L+ on spin 2.
%
% Example:
%
%    Sum_Lz=operator(spin_system,'Lz','13C');
%
% would return the sum of [Lz, ] commutation superoperators on all carbons in
% the system.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function A=operator(spin_system,operators,spins,operator_type)

% Preallocate the superoperator
A=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,spin_system.bas.nstates*spin_system.comp.nspins);

% Set the default type
if nargin==3
    operator_type='comm';
end

% Determine the output type
if ischar(spins)
    
    % If spins are supplied by name, parse the specification
    switch spins
        case 'all'
            spin_numbers=1:spin_system.comp.nspins;
        otherwise
            spin_numbers=find(strcmp(spin_system.comp.isotopes,spins));
    end
    
    % And return a sum of single-spin operators
    for n=1:numel(spin_numbers)
        A=A+operator(spin_system,operators,spin_numbers(n),operator_type);
    end
    
else

    % If spins are supplied by number, parse the specification
    [opspec,coeff]=human2opspec(spin_system,operators,spins);
    
    % And return the multi-spin order
    switch operator_type
        case 'left'
            A=coeff*p_superop(spin_system,opspec,'left');
        case 'right'
            A=coeff*p_superop(spin_system,opspec,'right');
        case 'comm'
            A=coeff*c_superop(spin_system,opspec);
        otherwise
            error('operator: unknown operator type.');
    end
    
end

end

% It is nice to know that the computer understands the problem. But I would
% like to understand it too.
%
% Eugene Wigner 
