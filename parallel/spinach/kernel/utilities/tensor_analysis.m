% Returns diagnostic information about an interacton tensor.
%
% ilya.kuprov@oerc.ox.ac.uk

function [eigvals,dcm,iso]=tensor_analysis(spin_system,tensor)

% Get the eigensystem
[dcm,eigvals]=eig(symmetrize(spin_system,full(tensor)));

% Sort the eigensystem
[~,index]=sort(abs(diag(eigvals)),1,'ascend');

% Rearrange the eigenvalues
eigvals=diag(eigvals);
eigvals=eigvals(index);

% Rearrange the eigenvectors
dcm=dcm(:,index);

% Kill the inversion component
dcm=dcm*det(dcm);

% Prefer upper half-space for positive directions
if dcm(3,3)<0
    dcm(:,3)=-dcm(:,3); dcm(:,1)=-dcm(:,1);
end

% Compute the isotropic part
iso=mean(eigvals);

end

% "He who dares not offend cannot be honest." - Thomas Paine

