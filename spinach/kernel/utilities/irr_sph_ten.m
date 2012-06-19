% Returns a cell array of single-spin irreducible spherical tensor operators T(k,m).
% 
% Two call formats are currently available. A two-argument call
%
%                               T=irr_sph_ten(mult,k)
%
% where 'mult' is the multiplicity of the spin in question and 'k' is the irreducible
% spherical tensor rank required, returns a cell array of tensors of that rank in the
% order of decreasing projection. A single argument call
%
%                               T=irr_sph_ten(mult)
%
% Produces tensors of all ranks and concatenates them in the order of ascending rank.
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function T=irr_sph_ten(mult,k)

switch nargin
    
    case 1
        
        % Generate tensors of all possible ranks and concatenate them into a single cell array
        T=cell(mult^2,1);
        for n=0:(mult-1)
            T((n^2+1):((n+1)^2))=irr_sph_ten(mult,n);
        end
        
    case 2
        
        % Generate tensors of the rank specified by the user
        if (k<0)||(k>(mult-1))
            % For non-physical ranks, bomb out
            error('irr_sph_ten: spherical tensor rank is outside the physical range.');
        elseif k==0
            % For zero rank, return the unit matrix
            T={speye(mult)};
        else
            
            % Get the Pauli matrices
            L=pauli(mult);
            
            % Preallocate the cell array
            T=cell(2*k+1,1);
            
            % Get the top state
            T{1}=(((-1)^k)*sqrt(factorial(2*k))/((2^k)*factorial(k)))*L.p^k;
            
            % Apply sequential lowering using Racah's commutation rule
            for n=2:(2*k+1)
                q=k-n+2;
                T{n}=(1/sqrt((k+q)*(k-q+1)))*(L.m*T{n-1}-T{n-1}*L.m);
            end
            
        end
        
end

end

% I swear by my life and my love of it that I will never live for the sake
% of another man, nor ask another man to live for mine.
%
% Ayn Rand, "Atlas Shrugged"