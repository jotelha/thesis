% Converts linear indexing state specification to L,M indexing. In
% the linear indexing convention, the states are listed in the order
% of increasing L rank, and, within ranks, in the order of decreas-
% ing M projection.
%
% WARNING - zero base indexing, that is: 
%             
%                           I=0 -> (L=0,M=0)
%                           I=1 -> (L=1,M=1)
%                           I=2 -> (L=1,M=0), et cetera...
% 
% Arrays of any dimension are accepted as arguments.
%
% ilya.kuprov@oerc.ox.ac.uk

function [L,M]=lin2lm(I)

% Make sure the input is valid
if any(nonzeros(I)<0); error('lin2lm: invalid linear index.'); end

% Get the ranks and projections
L=fix(sqrt(I)); M=L.^2+L-I;

% Make sure the conversion is correct
if nnz(lm2lin(L,M)~=I)>0
    error('lin2lm: IEEE arithmetic breakdown, please contact the developer.');
end

end

% Arrogance on the part of the meritorious is even more offensive
% to us than the arrogance of those without merit: for merit itself
% is offensive.
%
% Friedrich Nietzsche

