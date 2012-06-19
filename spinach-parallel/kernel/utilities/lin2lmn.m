% Converts linear indexing specification of a Wigner function to L,M,N
% indexing. In the linear indexing convention, the Wigner functions are
% listed in the order of increasing L rank. Within each L rank, the func-
% tions are listed in the order of decreasing left index, and, for each
% left index, in the order of decreasing right index.
%
% Wigner functions are enumerated using base one indexing, that is: 
%             
%                       I=1 -> (L=0,M=0,N=0)
%                       I=2 -> (L=1,M=1,N=1)
%                       I=3 -> (L=1,M=1,N=0), et cetera...
% 
% Arrays of any dimension are accepted as arguments.
%
% ilya.kuprov@oerc.ox.ac.uk

function [L,M,N]=lin2lmn(I)

% Make sure the input is valid
if any(nonzeros(I)<1); error('lin2lmn: invalid linear index.'); end

% Get the rank
big_root=(27*I+sqrt(729*I.^2-3)).^(1/3);
L=ceil((3^(1/3)+big_root.^2)./(2*(3^(2/3))*big_root)-1);

% Get the left index
rank_page_position=I-(4*L.^3-L)/3-1;
M=L-fix(rank_page_position./(2*L+1));

% Get the right index
N=L+(2*L+1).*(L-M)-rank_page_position;

% Make sure the conversion is correct
if nnz(lmn2lin(L,M,N)~=I)>0
    error('lin2lmn: IEEE arithmetic breakdown, please contact the developer.');
end

end

% "Who gave you the authority to decide which colour the Spinach logo
% was going to be?"
%
% Kelly-Anne Ferguson to IK, in April 2011.

