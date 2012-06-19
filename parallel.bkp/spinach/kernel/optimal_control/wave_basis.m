% Common basis sets for the expansion of pulse waveforms. Returns the 
% waveform basis functions as columns of a matrix.
%
% Because the resulting waveforms are discretized, they are not preci-
% sely orthogonal under the standard scalar multiplication. An extra
% orthogonalization step is therefore applied to make them orthogonal
% as vectors.
%
% ilya.kuprov@oerc.ox.ac.uk

function basis_waves=wave_basis(basis_type,n_functions,n_steps)

switch basis_type
    
    case 'sine_waves'
        
        % Preallocate the array
        basis_waves=zeros(n_functions,n_steps);
        
        % Fill the array
        for n=1:n_functions
            basis_waves(n,:)=sin(n*linspace(0,pi,n_steps));
        end
        
        % Orthogonalize the array
        basis_waves=orth(basis_waves');
        
    case 'cosine_waves'
        
        % Preallocate the array
        basis_waves=zeros(n_functions,n_steps);
        
        % Fill the array
        for n=1:n_functions
            basis_waves(n,:)=cos((n-1)*linspace(0,pi,n_steps));
        end
        
        % Orthogonalize the array
        basis_waves=orth(basis_waves');
        
    case 'all_waves'
        
        % Preallocate the array
        basis_waves=zeros(2*n_functions,n_steps);
        
        % Fill the array
        for n=1:n_functions
            basis_waves(2*n-1,:)=cos((n-1)*linspace(-pi,pi,n_steps));
            basis_waves(2*n,:)=sin(n*linspace(-pi,pi,n_steps));
        end
        
        % Orthogonalize the array
        basis_waves=orth(basis_waves');
        
    case 'legendre'
        
        % Preallocate the array
        basis_waves=zeros(n_functions,n_steps);
        
        % Fill the array
        for n=1:n_functions
            basis_waves(n,:)=mfun('P',n-1,linspace(-1,1,n_steps));
            basis_waves(n,:)=basis_waves(n,:)/norm(basis_waves(n,:));
        end
        
        % Orthogonalize the array
        basis_waves=orth(basis_waves');
        
end  
        
end

% Self-esteem is the reputation we acquire with ourselves. 
%
% Nathaniel Branden

