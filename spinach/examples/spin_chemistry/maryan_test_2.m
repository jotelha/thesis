% A template simulation for David Wilcox.
%
% ilya.kuprov@oerc.ox.ac.uk

function maryan_test_2()

% Set the isotopes
sys.isotopes={'E','E','14N','14N','1H'};

% Set the core basis set
bas.mode='IK-0';
bas.level=3;

% Set trajectory level restriction
% sys.tols.zte_tol=1e-3;

% Set the rotation matrices
R1=[0.438 0.8655 -0.2432; 0.8981 -0.4097 0.1595;-0.0384 0.2883 0.9568];
R2=[0.9703 -0.2207 0.0992; 0.2383 0.9426 -0.234;-0.0419 0.2506 0.9672];
R3=[0.9819 0.1883 -0.0203; -0.0348 0.285 0.9579; -0.1861 0.9398 -0.2864];

% Set the eigenvalues (Gauss)
A1=[-1.049  0 0; 0 -0.996 0; 0 0 13.826];
A2=[-0.305  0 0; 0 -0.222 0; 0 0 6.872];
A3=[-13.850 0 0; 0 -9.372 0; 0 0 0.143];

% Set the coupling tensors
inter.coupling.matrix=cell(5);
inter.coupling.matrix{1,3}=1e6*gauss2mhz(R1'*A1*R1);
inter.coupling.matrix{2,4}=1e6*gauss2mhz(R2'*A2*R2);
inter.coupling.matrix{1,5}=1e6*gauss2mhz(R3'*A3*R3);

% Set the zeeman interactions
inter.zeeman.scalar={2.0023 2.0023 0 0 0};

% Run Spinach
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
M=maryan(spin_system,50e-6,50e6,'icosahedral_642_vertices');

% Plot the answer
tri=load([spin_system.sys.root_dir spin_system.sys.slash 'exp' ...
                                   spin_system.sys.slash 'grids' ...
                                   spin_system.sys.slash 'icosahedral_642_triangles.txt'],'ASCII');
[X,Y,Z]=sph2cart(M(:,2),pi/2-M(:,1),M(:,3)-mean(M(:,3)));
trisurf(tri+1,X,Y,Z,M(:,3)-mean(M(:,3)));

end