% Tolerances and fundamental constants.
%
% ilya.kuprov@oerc.ox.ac.uk
% hannah.hogben@chem.ox.ac.uk

function spin_system=tolerances(spin_system,sys)

% Liouvillian zero tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'liouv_zero')
    spin_system.tols.liouv_zero=sys.tols.liouv_zero;
    report(spin_system,[pad('tolerances: Liouvillian zero tolerance',80) pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.liouv_zero=1e-7;
    report(spin_system,[pad('tolerances: Liouvillian zero tolerance',80) pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (safe default)']);
end

% Relaxation superoperator zero tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'rlx_zero')
    spin_system.tols.rlx_zero=sys.tols.rlx_zero;
    report(spin_system,[pad('tolerances: relaxation superoperator zero tolerance',80) pad(num2str(spin_system.tols.rlx_zero,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.rlx_zero=1e-7;
    report(spin_system,[pad('tolerances: relaxation superoperator zero tolerance',80) pad(num2str(spin_system.tols.rlx_zero,'%0.8e'),20) ' (safe default)']);
end

% Sparse matrix exponentiation method
if isfield(sys,'tols')&&isfield(sys.tols,'exponentiation')
    spin_system.tols.exponentiation=sys.tols.exponentiation;
    report(spin_system,[pad('tolerances: sparse matrix exponentiation method',80) pad(spin_system.tols.exponentiation,20) ' (user-specified)']);
else
    spin_system.tols.exponentiation='taylor';
    report(spin_system,[pad('tolerances: sparse matrix exponentiation method',80) pad(spin_system.tols.exponentiation,20) ' (safe default)']);
end

% Final zero tolerance in the exponential propagator
if isfield(sys,'tols')&&isfield(sys.tols,'prop_zero')
    spin_system.tols.prop_zero=sys.tols.prop_zero;
    report(spin_system,[pad('tolerances: exponential propagator zero tolerance',80) pad(num2str(spin_system.tols.prop_zero,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.prop_zero=1e-8;
    report(spin_system,[pad('tolerances: exponential propagator zero tolerance',80) pad(num2str(spin_system.tols.prop_zero,'%0.8e'),20) ' (safe default)']);
end

% Norm tolerance for the series terms in the exponential propagator
if isfield(sys,'tols')&&isfield(sys.tols,'prop_norm')
    spin_system.tols.prop_norm=sys.tols.prop_norm;
    report(spin_system,[pad('tolerances: norm tolerance for series terms of exponential propagator',80) pad(num2str(spin_system.tols.prop_norm,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.prop_norm=1e-9;
    report(spin_system,[pad('tolerances: norm tolerance for series terms of exponential propagator',80) pad(num2str(spin_system.tols.prop_norm,'%0.8e'),20) ' (safe default)']);
end

% Zero tolerance for the series terms in the exponential propagator
if isfield(sys,'tols')&&isfield(sys.tols,'prop_chop')
    spin_system.tols.prop_chop=sys.tols.prop_chop;
    report(spin_system,[pad('tolerances: zero tolerance for series terms of exponential propagator',80) pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.prop_chop=1e-10;
    report(spin_system,[pad('tolerances: zero tolerance for series terms of exponential propagator',80) pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (safe default)']);
end

% Final zero tolerance in the derivative propagator
if isfield(sys,'tols')&&isfield(sys.tols,'derivative_prop_zero')
    spin_system.tols.derivative_prop_zero=sys.tols.derivative_prop_zero;
    report(spin_system,[pad('tolerances: final zero tolerance in the derivative propagator',80) pad(num2str(spin_system.tols.derivative_prop_zero,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.derivative_prop_zero=1e-8;
    report(spin_system,[pad('tolerances: final zero tolerance in the derivative propagator',80) pad(num2str(spin_system.tols.derivative_prop_zero,'%0.8e'),20) ' (safe default)']);
end

% Norm tolerance for the series terms in the derivative propagator
if isfield(sys,'tols')&&isfield(sys.tols,'derivative_prop_norm')
    spin_system.tols.derivative_prop_norm=sys.tols.derivative_prop_norm;
    report(spin_system,[pad('tolerances: norm tolerance for the series terms of the derivative propagator',80) pad(num2str(spin_system.tols.derivative_prop_norm,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.derivative_prop_norm=1e-9;
    report(spin_system,[pad('tolerances: norm tolerance for the series terms of the derivative propagator',80) pad(num2str(spin_system.tols.derivative_prop_norm,'%0.8e'),20) ' (safe default)']);
end

% Zero tolerance for the series terms in the derivative propagator
if isfield(sys,'tols')&&isfield(sys.tols,'derivative_prop_chop')
    spin_system.tols.derivative_prop_chop=sys.tols.derivative_prop_chop;
    report(spin_system,[pad('tolerances: zero tolerance for the series terms of the derivative propagator',80) pad(num2str(spin_system.tols.derivative_prop_chop,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.derivative_prop_chop=1e-10;
    report(spin_system,[pad('tolerances: zero tolerance for the series terms of the derivative propagator',80) pad(num2str(spin_system.tols.derivative_prop_chop,'%0.8e'),20) ' (safe default)']);
end

% Krylov subspace dimension for the Krylov propagator
if isfield(sys,'tols')&&isfield(sys.tols,'krylov_dim')
    spin_system.tols.krylov_dim=sys.tols.krylov_dim;
    report(spin_system,[pad('tolerances: subspace dimension for Krylov propagator',80) pad(num2str(spin_system.tols.krylov_dim),20) ' (user-specified)']);
else
    spin_system.tols.krylov_dim=30;
    report(spin_system,[pad('tolerances: subspace dimension for Krylov propagator',80) pad(num2str(spin_system.tols.krylov_dim),20) ' (safe default)']);
end

% Krylov method switchover
if isfield(sys,'tols')&&isfield(sys.tols,'krylov_switchover')
    spin_system.tols.krylov_switchover=sys.tols.krylov_switchover;
    report(spin_system,[pad('tolerances: Krylov propagation to be used for nnz(L) above',80) pad(num2str(spin_system.tols.krylov_switchover),20) ' (user-specified)']);
else
    spin_system.tols.krylov_switchover=50000;
    report(spin_system,[pad('tolerances: Krylov propagation to be used for nnz(L) above',80) pad(num2str(spin_system.tols.krylov_switchover),20) ' (safe default)']);
end

% ZTE zero track tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'zte_tol')
    spin_system.tols.zte_tol=sys.tols.zte_tol;
    report(spin_system,[pad('tolerances: ZTE zero track tolerance',80) pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.zte_tol=1e-24;
    report(spin_system,[pad('tolerances: ZTE zero track tolerance',80) pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (safe default)']);
end

% ZTE sample length
if isfield(sys,'tols')&&isfield(sys.tols,'zte_nsteps')
    spin_system.tols.zte_nsteps=sys.tols.zte_nsteps;
    report(spin_system,[pad('tolerances: ZTE sample length',80) pad(num2str(spin_system.tols.zte_nsteps),20) ' (user-specified)']);
else
    spin_system.tols.zte_nsteps=16;
    report(spin_system,[pad('tolerances: ZTE sample length',80) pad(num2str(spin_system.tols.zte_nsteps),20) ' (safe default)']);
end

% ZTE state vector density threshold 
if isfield(sys,'tols')&&isfield(sys.tols,'zte_maxden')
    spin_system.tols.zte_maxden=sys.tols.zte_maxden;
    report(spin_system,[pad('tolerances: ZTE off for state vector densities in excess of',80) pad(num2str(spin_system.tols.zte_maxden),20) ' (user-specified)']);
else
    spin_system.tols.zte_maxden=0.5;
    report(spin_system,[pad('tolerances: ZTE off for state vector densities in excess of',80) pad(num2str(spin_system.tols.zte_maxden),20) ' (safe default)']);
end

% Proximity tolerance for the distance-based clusterizer
if isfield(sys,'tols')&&isfield(sys.tols,'prox_cutoff')
    spin_system.tols.prox_cutoff=sys.tols.prox_cutoff;
    report(spin_system,[pad('tolerances: proximity cut-off distance (Angstrom)',80) pad(num2str(spin_system.tols.prox_cutoff,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.prox_cutoff=4.0;
    report(spin_system,[pad('tolerances: proximity cut-off distance (Angstrom)',80) pad(num2str(spin_system.tols.prox_cutoff,'%0.8e'),20) ' (safe default)']);
end

% Interaction clean-up tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'inter_cutoff')
    spin_system.tols.inter_cutoff=sys.tols.inter_cutoff;
    report(spin_system,[pad('tolerances: interaction cutoff tolerance',80) pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.inter_cutoff=1e-5;
    report(spin_system,[pad('tolerances: interaction cutoff tolerance',80) pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (safe default)']);
end

% Basis printing hush tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'basis_hush')
    spin_system.tols.basis_hush=sys.tols.basis_hush;
    report(spin_system,[pad('tolerances: basis hush tolerance',80) pad(num2str(spin_system.tols.basis_hush),20) ' (user-specified)']);
else
    spin_system.tols.basis_hush=256;
    report(spin_system,[pad('tolerances: basis hush tolerance',80) pad(num2str(spin_system.tols.basis_hush),20) ' (safe default)']);
end

% Interaction tensor symmetry tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'inter_sym')
    spin_system.tols.inter_sym=sys.tols.inter_sym;
    report(spin_system,[pad('tolerances: interaction tensor symmetry tolerance',80) pad(num2str(spin_system.tols.inter_sym,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.inter_sym=1e-5;
    report(spin_system,[pad('tolerances: interaction tensor symmetry tolerance',80) pad(num2str(spin_system.tols.inter_sym,'%0.8e'),20) ' (safe default)']);
end

% Path tracing tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'path_trace')
    spin_system.tols.path_trace=sys.tols.path_trace;
    report(spin_system,[pad('tolerances: path tracing tolerance',80) pad(num2str(spin_system.tols.path_trace,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.path_trace=1e-7;
    report(spin_system,[pad('tolerances: path tracing tolerance',80) pad(num2str(spin_system.tols.path_trace,'%0.8e'),20) ' (safe default)']);
end

% Path drop tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'path_drop')
    spin_system.tols.path_drop=sys.tols.path_drop;
    report(spin_system,[pad('tolerances: path drop tolerance',80) pad(num2str(spin_system.tols.path_drop,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.path_drop=1e-10;
    report(spin_system,[pad('tolerances: path drop tolerance',80) pad(num2str(spin_system.tols.path_drop,'%0.8e'),20) ' (safe default)']);
end

% Irrep population tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'irrep_drop')
    spin_system.tols.irrep_drop=sys.tols.irrep_drop;
    report(spin_system,[pad('tolerances: irrep population tolerance',80) pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.irrep_drop=1e-10;
    report(spin_system,[pad('tolerances: irrep population tolerance',80) pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (safe default)']);
end

% Sparse algebra tolerance on dimension
if isfield(sys,'tols')&&isfield(sys.tols,'small_matrix')
    spin_system.tols.small_matrix=sys.tols.small_matrix;
    report(spin_system,[pad('tolerances: sparse algebra to be used for matrix dimension above',80) pad(num2str(spin_system.tols.small_matrix),20) ' (user-specified)']);
else
    spin_system.tols.small_matrix=200;
    report(spin_system,[pad('tolerances: sparse algebra to be used for matrix dimension above',80) pad(num2str(spin_system.tols.small_matrix),20) ' (safe default)']);
end

% Sparse algebra tolerance on density
if isfield(sys,'tols')&&isfield(sys.tols,'dense_matrix')
    spin_system.tols.dense_matrix=sys.tols.dense_matrix;
    report(spin_system,[pad('tolerances: sparse algebra to be used for matrix density below',80) pad(num2str(spin_system.tols.dense_matrix),20) ' (user-specified)']);
else
    spin_system.tols.dense_matrix=0.25;
    report(spin_system,[pad('tolerances: sparse algebra to be used for matrix density below',80) pad(num2str(spin_system.tols.dense_matrix),20) ' (safe default)']);
end

% Relative accuracy of the elements of Redfield superoperator
if isfield(sys,'tols')&&isfield(sys.tols,'rlx_integration')
    spin_system.tols.rlx_integration=sys.tols.rlx_integration;
    report(spin_system,[pad('tolerances: relative accuracy of the elements of Redfield superoperator',80) pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (user-specified)']);
else
    spin_system.tols.rlx_integration=1e-4;
    report(spin_system,[pad('tolerances: relative accuracy of the elements of Redfield superoperator',80) pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (safe default)']);
end

% Lebedev grid rank
if isfield(sys,'tols')&&isfield(sys.tols,'grid_rank')
    spin_system.tols.grid_rank=sys.tols.grid_rank;
    report(spin_system,[pad('tolerances: Lebedev grid rank',80) pad(num2str(spin_system.tols.grid_rank),20) ' (user-specified)']);
else
    spin_system.tols.grid_rank=131;
    report(spin_system,[pad('tolerances: Lebedev grid rank',80) pad(num2str(spin_system.tols.grid_rank),20) ' (safe default)']);
end

% Algorithm selection for propagator derivatives
if isfield(sys,'tols')&&isfield(sys.tols,'dP_method')
    spin_system.tols.dP_method=sys.tols.dP_method;
    report(spin_system,[pad('tolerances: propagator differentiation algorithm',80) pad(spin_system.tols.dP_method,20) ' (user-specified)']);
else
    spin_system.tols.dP_method='hausdorff';
    report(spin_system,[pad('tolerances: propagator differentiation algorithm',80) pad(spin_system.tols.dP_method,20) ' (safe default)']);
end

% Hausdorff series accuracy for propagator derivatives
if isfield(sys,'tols')&&isfield(sys.tols,'dP_order')
    spin_system.tols.dP_order=sys.tols.dP_order;
    report(spin_system,[pad('tolerances: Hausdorff series order for propagator gradients',80) pad(num2str(spin_system.tols.dP_order),20) ' (user-specified)']);
else
    spin_system.tols.dP_order=Inf;
    report(spin_system,[pad('tolerances: Hausdorff series order for propagator gradients',80) pad(num2str(spin_system.tols.dP_order),20) ' (safe default)']);
end

% Norm estimation tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'normest_tol')
    spin_system.tols.normest_tol=sys.tols.normest_tol;
    report(spin_system,[pad('tolerances: norm estimation tolerance',80) pad(num2str(spin_system.tols.normest_tol),20) ' (user-specified)']);
else
    spin_system.tols.normest_tol=5e-4;
    report(spin_system,[pad('tolerances: norm estimation tolerance',80) pad(num2str(spin_system.tols.normest_tol),20) ' (safe default)']);
end

% Fundamental constants
spin_system.tols.hbar=1.054571628e-34;
report(spin_system,[pad('tolerances: Planck constant (hbar)',80) pad(num2str(spin_system.tols.hbar,'%0.8e'),20)]);
spin_system.tols.kbol=1.3806503e-23;
report(spin_system,[pad('tolerances: Boltzmann constant (k)',80) pad(num2str(spin_system.tols.kbol,'%0.8e'),20)]);
spin_system.tols.freeg=2.0023193043622;
report(spin_system,[pad('tolerances: free electron g-factor',80) pad(num2str(spin_system.tols.freeg,'%0.8e'),20)]);

end

% "Man once surrendering his reason, has no remaining guard against
% absurdities the most monstrous, and like a ship without rudder, is the
% sport of every wind." - Thomas Jefferson

