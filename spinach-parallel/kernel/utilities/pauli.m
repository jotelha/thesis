% Pauli matrices for a spin of user-specified multiplicity.
%
% ilya.kuprov@oerc.ox.ac.uk

function pauli=pauli(mult)

% Make sure the input is valid
if mult<0||(ceil(mult)~=floor(mult))||(numel(mult)~=1)||(imag(mult)~=0)
    error('pauli: invalid spin multiplicity specification.');
end

% Get the Pauli matrices
if mult==2
    % Spin-half matrices are hard-coded for speed
    pauli.p=sparse([0 1; 0 0]);
    pauli.m=sparse([0 0; 1 0]);
    pauli.z=sparse([0.5 0; 0 -0.5]);
    pauli.x=0.5*(pauli.p+pauli.m); pauli.y=-0.5*1i*(pauli.p-pauli.m);
else
    % Everything else goes through the standard procedure
    spin=(mult-1)/2;
    prjs=((mult-1):-1:0)-spin;
    pauli.p=spdiags(sqrt(spin*(spin+1)-prjs.*(prjs+1))',1,mult,mult);
    pauli.m=spdiags(sqrt(spin*(spin+1)-prjs.*(prjs-1))',-1,mult,mult);
    pauli.x=0.5*(pauli.p+pauli.m);
    pauli.y=-0.5*1i*(pauli.p-pauli.m);
    pauli.z=spdiags(prjs',0,mult,mult);
end
    
end

% Any refusal to recognize reality, for any reason whatever, has
% distrastrous consequences. There are no evil thoughts except one - the
% refusal to think. Don't ignore your own desires... Don't sacrifice them.
% Examine their cause.
%
% Ayn Rand, "Atlas Shrugged"