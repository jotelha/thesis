% Anisotropic singlet-singlet MARY experiment with exponential
% recombination function.
%
% field:    magnetic induction in Tesla
% kinetics: singlet recombination rate constant in Hz
% grid:     spherical grid to be used
%
% The code runs in parallel over the orientations.
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function M=maryan(spin_system,field,kinetics,grid)

% Show the banner
banner(spin_system,'sequence_banner');

% Get the spherical grid
grid=load([spin_system.sys.root_dir spin_system.sys.slash 'exp' ...
                                    spin_system.sys.slash 'grids' ...
                                    spin_system.sys.slash grid '.txt'],'ASCII');
grid_size=size(grid,1);

% Get the fields into angular frequencies
gamma=spin('E'); field=field*gamma;

% Get the two-electron singlet state
electrons=find(strcmp(spin_system.comp.isotopes,'E'));
if length(electrons)~=2
    error('maryan: number of electrons not equal to two.');
end
S=singlet(spin_system,electrons(1),electrons(2));

% Generate the field operators
report(spin_system,'maryan: generating the field operators...');
E1z=operator(spin_system,'Lz',electrons(1));
E2z=operator(spin_system,'Lz',electrons(2));

% Get the rotation basis
report(spin_system,'maryan: building the rotation basis...');
[Iso,Q]=h_superop(secularity(spin_system,'lowfield'));

% Run the state space restriction
report(spin_system,'maryan: trying to reduce the problem dimension...');
projectors=reduce(spin_system,Iso+orientation(Q,[0 pi/11 pi/7])+field*(E1z+E2z),S);

% Determine the Euler angles
[azimuth,elevation]=cart2sph(grid(:,1),grid(:,2),grid(:,3));
theta=pi/2-elevation; phi=azimuth;

% Preallocate the result array
M=zeros(grid_size,1);

% Loop over the independent subspaces
for subs=1:length(projectors)
    
    % Project into the current subspace
    Iso_subs=projectors{subs}'*Iso*projectors{subs};
    Q_subs=cell(size(Q));
    for n=1:5
        for k=1:5
            Q_subs{n,k}=projectors{subs}'*Q{n,k}*projectors{subs};
        end
    end         
    E1z_subs=projectors{subs}'*E1z*projectors{subs};
    E2z_subs=projectors{subs}'*E2z*projectors{subs};
    S_subs=projectors{subs}'*S;

    % Loop over the grid points
    parfor n=1:grid_size
    
        % Get the Liouvillian for the current orientation
        L=Iso_subs+orientation(Q_subs,[phi(n) theta(n) 0])+field*(E1z_subs+E2z_subs);
        
        % Compute MARY
        M(n)=M(n)+kinetics*real(S_subs'*((1i*L+kinetics*speye(size(L)))\S_subs));
        
    end
    
end

% Return the result
M=[theta phi M];

end