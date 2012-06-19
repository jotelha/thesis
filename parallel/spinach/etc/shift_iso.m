% Replaces the isotropic parts of interaction tensors with user-supplied values.
% This is useful for correcting DFT calculations, where the anisotropy is usually
% satisfactory, but the isotropic art often is not. Arguments:
%
%      tensors      - a cell array of interaction tensors as 3x3 matrices
%      spin_numbers - a vector containing the numbers of spins in the tensors
%                     array that should have the isotropic values replaced
%      new_iso      - a vector containing the new isotropic parts
%
% luke.edwards@chem.ox.ac.uk
% ilya.kuprov@oerc.ox.ac.uk

function tensors=shift_iso(tensors,spin_numbers,new_iso)

% Loop over the tensors
for n=1:numel(spin_numbers)
    
    % Pull out the anisotropy
    [~,~,rank2]=mat2sphten(tensors{spin_numbers(n)});
    
    % Rebuild with the new isotropic part
    tensors{spin_numbers(n)}=sphten2mat([],[],rank2)+new_iso(n)*eye(3);
    
end

end

% http://news.bbc.co.uk
% http://www.truecrypt.org/downloads
% http://www.torproject.org/download/download-easy.html
% http://kpvz7ki2v5agwt35.onion

