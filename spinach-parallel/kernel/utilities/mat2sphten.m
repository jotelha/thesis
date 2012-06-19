% Converts a 3x3 interaction matrix into the irreducible spherical tensor
% notation: one rank 0 component, three rank 1 components and five rank 2
% components to the total of nine independent components.
%
% The components are listed in the following order:
%
% rank 0: (0,0)
% rank 1: (1,1) (1,0) (1,-1)
% rank 2: (2,2) (2,1) (2,0) (2,-1) (2,-2)
%
% See Table 1 in:
%       P.L. Nordio, U. Serge, J. Magn. Reson., 1977, 27, 465-473.
%
% ilya.kuprov@oerc.ox.ac.uk

function [rank0,rank1,rank2]=mat2sphten(M)

% Rank 0 component
rank0=-(1/sqrt(3))*trace(M);

% Rank 1 components
rank1(1)=-(1/2)*(M(3,1)-M(1,3)-1i*(M(3,2)-M(2,3)));
rank1(2)=-(1i/sqrt(2))*(M(1,2)-M(2,1));
rank1(3)=-(1/2)*(M(3,1)-M(1,3)+1i*(M(3,2)-M(2,3)));

% Rank 2 components
rank2(1)=+(1/2)*(M(1,1)-M(2,2)-1i*(M(1,2)+M(2,1)));
rank2(2)=-(1/2)*(M(1,3)+M(3,1)-1i*(M(2,3)+M(3,2)));
rank2(3)=(1/sqrt(6))*(2*M(3,3)-M(1,1)-M(2,2));
rank2(4)=+(1/2)*(M(1,3)+M(3,1)+1i*(M(2,3)+M(3,2)));
rank2(5)=+(1/2)*(M(1,1)-M(2,2)+1i*(M(1,2)+M(2,1)));

end

% "When I was 13 I think - Hewlett and Packard were my idols - I called
% up Bill Hewlett because he lived in Palo Alto and there were no unlisted
% numbers in the phonebook - which gives you a clue to my age. And he
% picked up the phone and I talked to him and I asked him if he'd give me
% some spare parts for something I was building called a frequency
% counter. And he did, but in addition to that he gave me something way
% more important. He gave me a job that summer - a summer job - at
% Hewlett Packard right here in Santa Clara off 280, in a division that
% built frequency counters. And I was in heaven."
%
% Steve Jobs

