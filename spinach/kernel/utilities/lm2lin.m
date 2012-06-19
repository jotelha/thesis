% Converts L,M spin state specification to linear indexing specification. In
% the linear indexing convention, the states are listed in the order of inc-
% reasing L rank, and, within ranks, in the order of decreasing M projection.
%
% WARNING - zero base indexing, that is: 
%             
%                           (L=0,M=0) -> I=0
%                           (L=1,M=1) -> I=1
%                           (L=1,M=0) -> I=2, et cetera...
%
% Arrays of any dimension are accepted as arguments.
%
% ilya.kuprov@oerc.ox.ac.uk

function I=lm2lin(L,M)

% Make sure the input is valid
if any(abs(M)>L); error('lm2lin: unacceptable projection number.'); end
if any(L<0); error('lm2lin: unacceptable total angular momentum.'); end
if any(size(L)~=size(M)); error('lm2lin: array dimensions are inconsistent.'); end

% Get the linear index
I=L.^2+L-M;

end

% A casual stroll through the lunatic asylum shows that faith does not
% prove anything.
%
% Friedrich Nietzsche

