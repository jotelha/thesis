% Converts the nine components of the irreducible spherical tensor
% representation of an interaction tensor into the Cartesian 
% representation with a 3x3 matrix. 
%
% Spherical tensor components should be listed in the following order:
%
% rank 0: (0,0)
% rank 1: (1,1) (1,0) (1,-1)
% rank 2: (2,2) (2,1) (2,0) (2,-1) (2,-2)
%
% See Table 1 in:
%       P.L. Nordio, U. Serge, J. Magn. Reson., 1977, 27, 465-473.
%
% ilya.kuprov@oerc.ox.ac.uk

function M=sphten2mat(rank0,rank1,rank2)

M=zeros(3);

% Rank 0 component
if ~isempty(rank0)
    M=M-(1/sqrt(3))*[1 0 0; 0 1 0; 0 0 1]*rank0;
end

% Rank 1 components
if exist('rank1','var')&&~isempty(rank1)
    M=M-(1/2)*[0 0 -1; 0 0 -1i; 1 1i 0]*rank1(1);
    M=M-(1/sqrt(8))*[0 -2i 0; 2i 0 0; 0 0 0]*rank1(2);
    M=M-(1/2)*[0 0 -1; 0 0 1i; 1 -1i 0]*rank1(3);
end

% Rank 2 components
if exist('rank2','var')&&~isempty(rank2)
    M=M+(1/2)*[1 1i 0; 1i -1 0; 0 0 0]*rank2(1);
    M=M-(1/2)*[0 0 1; 0 0 1i; 1 1i 0]*rank2(2);
    M=M+(1/sqrt(6))*[-1 0 0; 0 -1 0; 0 0 2]*rank2(3);
    M=M+(1/2)*[0 0 1; 0 0 -1i; 1 -1i 0]*rank2(4);
    M=M+(1/2)*[1 -1i 0; -1i -1 0; 0 0 0]*rank2(5);
end

end

% "When you're young, you look at television and think "there's a
% conspiracy - the networks have conspired to dumb us down". But when 
% you get a little older, you realize that's not true. The networks
% are in business to give people exactly what they want. That's a 
% far more depressing thought."
%
% Steve Jobs

