% Converts user-friendly descriptions of spin states and operators into the 
% formal description (opspec) used by Spinach kernel. Arguments:
%	
%  operators - a cell array of strings. Each string may be be 'E' (identity 
%              operator), 'Lz', 'L+', 'L-' (single spin orders) or 'Tl,m'
%              (higher spin orders as irreducible spherical tensors; l and
%              m are both integers).
%
%      spins - a cell array of integers specifying the numbers of spins on
%              which the operators given in the 'oper' argument operate.
%
% Example:
%
%   [opspec,coefficient]=human2opspec(spin_system,{'Lz','L+'},{1,2});
%
% would return the required opspec and coefficient to generate an LzL+ operator
% or state vector with Lz on the spins 1 and L+ on spin 2.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function [opspec,coefficient]=human2opspec(spin_system,operators,spins)

% Convert inputs to cell arrays
if ~iscell(operators), operators={operators}; end
if ~iscell(spins), spins={spins}; end

% Validate the input
if numel(operators)~=numel(spins)
    error('human2opspec: dimensions of `oper` and `spins` arguments do not match.')
end
if numel(unique(cell2mat(spins)))~=numel(spins)
    error('human2opspec: repetitions detected in `spins` input.')
end

% Start with an empty opspec and a unit coefficient
opspec=zeros(1,spin_system.comp.nspins);
coefficient=1;

% Parse operator selection.
for n=1:numel(operators)
    
    switch operators{n}
        case 'E'
            opspec(spins{n})=0;
            coefficient=1*coefficient;
        case 'L+'
            opspec(spins{n})=1;
            coefficient=-sqrt(2)*coefficient;
        case 'Lz'
            opspec(spins{n})=2;
            coefficient=1*coefficient;
        case 'L-'
            opspec(spins{n})=3;
            coefficient=sqrt(2)*coefficient;
        otherwise
            
            % Validate the irreducible spherical tensor input
            if isempty(regexp(operators{n},'^T([\+\-]?\d+),([\+\-]?\d+)$','once'))
                error('human2opspec: unrecognized operator specification.');
            end             
            
            % Extract the quantum numbers
            indices=textscan(operators{n},'T%n,%n');
            l=indices{1}; m=indices{2};
                        
            % Validate the quantum numbers
            if (l<0)||(abs(m)>l)
                error('human2opspec: invalid indices in irreducible spherical tensors.');
            end
            
            % Write the specification
            opspec(spins{n})=lm2lin(l,m);
            coefficient=1*coefficient;
            
    end
    
end

end

% God made the bulk; surfaces were invented by the Devil.
%
% Wolfgang Pauli



