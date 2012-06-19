% Builds the generators to be used in Hilbert space simulations.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function generator_array=hs_generators(spin_system)

% Generate the database record name
if ispc
    generator_file=[spin_system.sys.root_dir '\kernel\cache\generators_' strrep(num2str(spin_system.comp.mults),' ','_') '.mat'];
elseif isunix||ismac
    generator_file=[spin_system.sys.root_dir '/kernel/cache/generators_' strrep(num2str(spin_system.comp.mults),' ','_') '.mat'];
end

% Check if the record exists in the database
if exist(generator_file,'file')==2
    
    % Load generators from the database
    load(generator_file,'-mat');
    
else
    
    % Preallocate the cell array
    generator_array=struct('Lx',num2cell(zeros(1,spin_system.comp.nspins)),...
                           'Ly',num2cell(zeros(1,spin_system.comp.nspins)),...
                           'Lz',num2cell(zeros(1,spin_system.comp.nspins)),...
                           'Lp',num2cell(zeros(1,spin_system.comp.nspins)),...
                           'Lm',num2cell(zeros(1,spin_system.comp.nspins)));
    
    % Loop over spins in the system
    for n=1:spin_system.comp.nspins
        
        % Get Pauli matrices for the current spin
        L=pauli(spin_system.comp.mults(n));
        
        % Build the generators
        A=speye(prod(spin_system.comp.mults(1:(n-1))));
        B=speye(prod(spin_system.comp.mults((n+1):end)));
        generator_array(n).Lx=kron(kron(A,L.x),B);
        generator_array(n).Ly=kron(kron(A,L.y),B);
        generator_array(n).Lz=kron(kron(A,L.z),B);
        generator_array(n).Lp=generator_array(n).Lx+1i*generator_array(n).Ly;
        generator_array(n).Lm=generator_array(n).Lx-1i*generator_array(n).Ly;
        
    end
    
    % Store the generators for future reference
    save(generator_file,'generator_array');
    
    % Inform the user
    report(spin_system,'generators: generator cache updated.');
    
end

end

% "There used to be a simple story about Russian literature, that we
% thought the good writers were the ones who opposed the regime. Once we
% don't have that story about Russia as a competitor, or an enemy, it was
% much less clear to us what we should be interested in."
%
% Edwin Frank, the editor of NYRB Classics

