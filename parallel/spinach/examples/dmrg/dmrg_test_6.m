% Test of translation functions between Spinach representation of
% superoperators and MPOs. Illustrates the fact that the operators
% are very compressible, but the propagators are not.
%
% ilya.kuprov@oerc.ox.ac.uk

function dmrg_test_6()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','1H','1H','1H','1H'};

% Couplings
inter.coupling.scalar=cell(5);
for n=1:4
    inter.coupling.scalar{n,n+1}=10+n;
end

% Basis set
bas.mode='complete';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Test the Hamiltonian translation and compression
H=h_superop(secularity(spin_system,'nmr'));
cito=spinach2cito(spin_system,full(H));
mpo=cito2mpo(spin_system,cito,1e-12);
cito=mpo2cito(mpo_compress(spin_system,mpo,1e-6));
norm(H-cito2spinach(spin_system,cito),'inf')

% Test the propagator translation
P=propagator(spin_system,H,stepsize(H));
cito=spinach2cito(spin_system,full(P));
mpo=cito2mpo(spin_system,cito,1e-12);
cito=mpo2cito(mpo_compress(spin_system,mpo,1e-6));
norm(P-cito2spinach(spin_system,cito),'inf')

end

