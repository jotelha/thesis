% Converts a direction cosine matrix into second-rank Wigner function matrix.
%
% Rows and columns of the matrix are sorted by descending rank, that is
%    
%                         [D( 2,2)  ...  D( 2,-2)
%                            ...    ...    ...  
%                          D(-2,2)  ...  D(-2,-2)]
%
% The resulting Wigner matrix is to be used as v=W*v, where v is a column
% vector of irreducible components, listed vertically in the following
% order: (2,2) (2,1) (2,0) (2,-1) (2,-2)
%
% ilya.kuprov@oerc.ox.ac.uk

function wigner=dcm2wigner(dcm)

% Set the deviation tolerance
tol=1e-6;

% Make sure the DCM is unitary
if norm(dcm'*dcm-eye(3))>tol
    error('dcm2wigner: the DCM is not orthogonal.');
end

% Make sure the DCM has a unit determinant
if abs(det(dcm)-1)>tol
    error('dcm2wigner: the DCM is not a proper rotation.');
end

% Compute the A and B coefficients
A=sqrt(0.5*(dcm(1,1)+1i*dcm(1,2)-1i*dcm(2,1)+dcm(2,2)));
B=sqrt(0.5*(-dcm(1,1)+1i*dcm(1,2)+1i*dcm(2,1)+dcm(2,2)));

% Verify the amplitudes
err=abs(A*A'-0.5*(1+dcm(3,3)))+abs(B*B'-0.5*(1-dcm(3,3)));
if err>tol
    error('dcm2wigner: DCM does not pass the self-consistency check on amplitudes.');
end

% Verify the phases
err=abs(A*B+0.5*(dcm(1,3)-1i*dcm(2,3)))+abs(A*B'-0.5*(dcm(3,1)+1i*dcm(3,2)));
if err>tol
    A=-A;
end
err=abs(A*B+0.5*(dcm(1,3)-1i*dcm(2,3)))+abs(A*B'-0.5*(dcm(3,1)+1i*dcm(3,2)));
if err>tol
    error('dcm2wigner: DCM does not pass the self-consistency check on phases.');
end

% Compute the Wigner matrix
Z=A*A'-B*B';
wigner=[ A'^4                2*A'^3*B'          sqrt(6)*A'^2*B'^2      2*A'*B'^3           B'^4             ;
        -2*A'^3*B            A'^2*(2*Z-1)       sqrt(6)*A'*B'*Z        B'^2*(2*Z+1)        2*A*B'^3         ;
         sqrt(6)*A'^2*B^2   -sqrt(6)*A'*B*Z     0.5*(3*Z^2-1)          sqrt(6)*A*(B')*Z    sqrt(6)*A^2*B'^2 ;
        -2*A'*B^3            B^2*(2*Z+1)       -sqrt(6)*A*B*Z          A^2*(2*Z-1)         2*A^3*B'         ;
         B^4                -2*A*B^3            sqrt(6)*A^2*B^2       -2*A^3*B             A^4             ];

end

% The hallmark of a second rater is resentment for another man's achievement.
%
% Ayn Rand, "Atlas Shrugged"

