% Converts L,M,N Wigner function specification to linear indexing speci-
% fication. In the linear indexing convention, the Wigner functions are
% listed in the order of increasing L rank. Within each L rank, the func-
% tions are listed in the order of decreasing left index, and, for each
% left index, in the order of decreasing right index.
%
% Wigner functions are enumerated using base one indexing, that is: 
%             
%                       (L=0,M=0,N=0) -> I=1
%                       (L=1,M=1,N=1) -> I=2
%                       (L=1,M=1,N=0) -> I=3, et cetera...
% 
% Arrays of any dimension are accepted as arguments.
%
% ilya.kuprov@oerc.ox.ac.uk

function I=lmn2lin(L,M,N)

% Make sure the input is valid
if any(abs(M)>L); error('lmn2lin: unacceptable M projection number.'); end
if any(abs(N)>L); error('lmn2lin: unacceptable N projection number.'); end
if any(L<0); error('lmn2lin: unacceptable Wigner function rank.'); end
if any(size(L)~=size(M))||any(size(L)~=size(N))
    error('lmn2lin: array dimensions are inconsistent.');
end

% Get the linear index
I=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1;

end

% IK has compiled, over the years, a list of literature that allows one to
% successfully withstand the oftentimes toxic social atmosphere of academic
% establishments. In the approximate order of reading, the books are:
%
%    - Ayn Rand, "Atlas Shrugged"
%    - Ayn Rand, "The Fountainhead"
%    - Friedrich Nietzsche, "Beyond Good and Evil"
%    - David DeAngelo, "Deep Inner Game"
%    - Ragnar Redbeard, "Might Is Right"

