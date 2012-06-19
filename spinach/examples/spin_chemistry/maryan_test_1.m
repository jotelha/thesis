% Figure 7.10 from the DPhil thesis by Chris Rodgers.
%
% ilya.kuprov@oerc.ox.ac.uk

function maryan_test_1()

% Set the common parameters
sys.isotopes={'E','E','1H'};
bas.mode='complete';

% Set the couplings
inter.zeeman.scalar={2.0023 2.0023 0};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,3}=1e6*gauss2mhz([5 0 0; 0 4 0; 0 0 10]);

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
M=maryan(spin_system,50e-6,2e5,'icosahedral_10242_vertices');

% Plot the answer
tri=load([spin_system.sys.root_dir spin_system.sys.slash 'exp' ...
                                   spin_system.sys.slash 'grids' ...
                                   spin_system.sys.slash 'icosahedral_10242_triangles.txt'],'ASCII');
[X,Y,Z]=sph2cart(M(:,2),pi/2-M(:,1),M(:,3)-mean(M(:,3)));
trisurf(tri+1,X,Y,Z,M(:,3)-mean(M(:,3)));

end