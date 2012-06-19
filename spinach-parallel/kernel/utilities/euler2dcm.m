% Converts Euler angles to DCM.
%
% The resulting rotation matrix is to be used as follows:
%
%     v=R*v    (for 3x1 vectors)
%     A=R*A*R' (for 3x3 interaction tensors)
%
% ilya.kuprov@oerc.ox.ac.uk

function R=euler2dcm(arg1,arg2,arg3)

% Adapt to the input style
if nargin==1
    alpha=arg1(1); beta=arg1(2); gamma=arg1(3);
elseif nargin==3
    alpha=arg1; beta=arg2; gamma=arg3;
else
    error('euler2dcm: incorrect number of input arguments.');
end

% Build the individual rotation matrices
R_alpha=[ cos(alpha)  sin(alpha) 0;
         -sin(alpha)  cos(alpha) 0;
          0           0          1];
R_beta= [cos(beta)   0         -sin(beta);
         0           1          0;
         sin(beta)   0          cos(beta)];
R_gamma=[ cos(gamma)  sin(gamma) 0;
         -sin(gamma)  cos(gamma) 0;
          0           0          1];
     
% Build the directional cosine matrix
R=R_gamma*R_beta*R_alpha;

end

% I am so hungry for any sight of anyone who's able to do whatever it is
% he's doing!
%
% Ayn Rand, "Atlas Shrugged"

