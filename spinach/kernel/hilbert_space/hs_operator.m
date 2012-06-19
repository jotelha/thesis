% Generates Hilbert space product operators from their user-friendly descriptions.
% Arguments:
%	
%     states - a cell array of strings. Each string may be 'Lz', 'Lx' or 'Ly',
%              the states of all spins not mentioned explicitly are assumed
%              to be identity.
%
%      spins - a cell array of integers specifying the numbers of spins
%              to be generated in states given in the 'states' argument.
%              
% Example:
%
%    LzSy=hs_operator(spin_system,{'Lz','Ly'},{1,2});
%
% would return the LzSy operator with Lz on spin 1 and Ly on spin 2.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function A=hs_operator(spin_system,operators,spins)

% Determine the output type...
if ischar(spins)
    
    % Preallocate the superoperator...
    A_dim=prod(spin_system.comp.mults);
    A=sparse(A_dim,A_dim);

    % If spins are supplied by name, parse the specification...
    switch spins
        case 'all'
            spin_numbers=1:spin_system.comp.nspins;
        otherwise
            spin_numbers=find(strcmp(spin_system.comp.isotopes,spins));
    end
    
    % And return a sum of single-spin operators.
    for n=1:numel(spin_numbers)
        A=A+hs_operator(spin_system,operators,spin_numbers(n));
    end
    
else

    % If spins are supplied by number, parse the specification.
    [opspec,coeff]=human2opspec(spin_system,operators,spins);
    
    % Take the Kronecker product of all the states in the opspec to form
    % the multi-spin operators.
    A=1;
    k=1;
    unit_dim=1;
    while k<=size(opspec,2)
        if opspec(k)==0
            % Corresponds to T_{0,0}, i.e. a unit matrix
            unit_dim=unit_dim*spin_system.comp.mults(k);
        elseif unit_dim>1
            % kron(I_n,I_m)=I_{nm}
            [L,M]=lin2lm(opspec(k));
            ist=irr_sph_ten(spin_system.comp.mults(k),L);
            A=kron(A,kron(speye(unit_dim),ist{L-M+1}));
            unit_dim=1;
        else
            % Normal operation
            [L,M]=lin2lm(opspec(k));
            ist=irr_sph_ten(spin_system.comp.mults(k),L);
            A=kron(A,ist{L-M+1});
        end
        k=k+1;
    end
    
	% Multiply by the appropriate coefficient and account for when opspec 
	% ends with a string of zeroes.
    A=coeff*kron(A,speye(unit_dim));

end

end

% "...you see, God - whatever anyone chooses to call God - is one's highest
% conception of the highest possible. And whoever places his highest
% conception above his own possibility thinks very little of himself and
% his life. It's a rare gift, you know, to feel reverence for your own life
% and to want the best, the greatest, the highest possible, here, now, for
% your very own. To imagine a heaven and then not to dream of it, but to
% demand it."
%
% Ayn Rand, "We the Living"

