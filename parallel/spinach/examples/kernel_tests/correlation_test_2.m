% Computes rotational correlation functions using a Monte-Carlo method and
% compares them to the analytical results supplied by the Spinach kernel.
%
% The problem has been scaled into a unit time step. Typical call:
%
%    correlation_test_2(0.1,0.2,2,5,2,5)
%
% ilya.kuprov@oerc.ox.ac.uk

function correlation_test_2(sigma_ax,sigma_eq,k,m,p,q)

% Set the number of points
npoints=1e5;

% Start with a unit DCM
DCM=eye(3);

% Generate the angle track
angles=randn(3,npoints);

% Preallocate the Wigner function array
W=zeros(5,5,npoints);

% Loop over the Monte-Carlo steps
for n=1:npoints
    
    % Generate a random rotation
    R_gen=sigma_ax*[ 0  1  0; -1  0  0;  0  0  0]*angles(1,n)+...
          sigma_eq*[ 0  0  1;  0  0  0; -1  0  0]*angles(2,n)+...
          sigma_eq*[ 0  0  0;  0  0  1;  0 -1  0]*angles(3,n);
    R=expm(R_gen);
   
    % Add the rotation to the current one
    DCM=R*DCM;
    
    % Convert the current orientation into Wigner functions
    W(:,:,n)=dcm2wigner(DCM);
       
end

% Get the Monte-Carlo correlation function
cf_mc=(1/5)*ifftshift(real(crosscorr(squeeze(W(k,m,:)),squeeze(W(p,q,:)),100)));
plot(cf_mc(1:100),'r.'); hold on;

% Get the analytical correlation function
spin_system.rlx.tau_c=1./(3*[sigma_ax sigma_eq].^2);
cf_an=correlation_function(spin_system,k,m,p,q,0:99);
plot(cf_an,'b-');

end
   
   
    
