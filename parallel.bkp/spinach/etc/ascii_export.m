% Writes ASCII files with sparse representations of the user-specified operators.
%
% ilya.kuprov@oerc.ox.ac.uk

function ascii_export(spin_system,operators,states)

% Dump the isotropic Hamiltonian commutation superoperator
if any(strcmp(operators,'H'))
    H=h_superop(spin_system);
    [i,j,v]=find(H); A=[i j real(v)]; B=[i j imag(v)];
    save('h_iso_real.txt','A','-ASCII');
    save('h_iso_imag.txt','B','-ASCII');
end

% Dump the rotational decomposition of the anisotropic part of the Hamiltonian commutatioin superoperator
if any(strcmp(operators,'Q'));
    [H,Q]=h_superop(spin_system);
    for n=1:5
        for k=1:5
            [i,j,v]=find(Q{n,k}); A=[i j real(v)]; B=[i j imag(v)];
            save(['h_aniso_' num2str(n) '_' num2str(k) '_real.txt'],'A','-ASCII');
            save(['h_aniso_' num2str(n) '_' num2str(k) '_imag.txt'],'B','-ASCII');
        end
    end
end

% Dump the relaxation superoperator
if any(strcmp(operators,'R'))
    R=r_superop(spin_system);
    [i,j,v]=find(R); A=[i j real(v)]; B=[i j imag(v)];
    save('rlx_real.txt','A','-ASCII');
    save('rlx_imag.txt','B','-ASCII');
end

% Dump the control superoperators
if any(strcmp(operators,'controls'));
    for n=1:spin_system.nspins
        Lm=sqrt(2)*c_superop(spin_system,sparse(1,n,3,1,spin_system.nspins,1));
        [i,j,v]=find(Lm); A=[i j real(v)]; B=[i j imag(v)];
        save(['Lm_superop_spin_' num2str(n) '_real.txt'],'A','-ASCII');
        save(['Lm_superop_spin_' num2str(n) '_imag.txt'],'B','-ASCII');
        [i,j,v]=find(Lm'); A=[i j real(v)]; B=[i j imag(v)];
        save(['Lp_superop_spin_' num2str(n) '_real.txt'],'A','-ASCII');
        save(['Lp_superop_spin_' num2str(n) '_imag.txt'],'B','-ASCII');
        [i,j,v]=find((Lm'*Lm-Lm*Lm')/2); A=[i j real(v)]; B=[i j imag(v)];
        save(['Lz_superop_spin_' num2str(n) '_real.txt'],'A','-ASCII');
        save(['Lz_superop_spin_' num2str(n) '_imag.txt'],'B','-ASCII');
    end
end

% Dump the states
if any(strcmp(states,'essential'));
    for n=1:spin_system.nspins
        Lm=sqrt(2)*statevec(spin_system,sparse(1,n,3,1,spin_system.nspins,1));
        [i,j,v]=find(Lm); A=[i j real(v)]; B=[i j imag(v)];
        save(['Lm_state_spin_' num2str(n) '_real.txt'],'A','-ASCII');
        save(['Lm_state_spin_' num2str(n) '_imag.txt'],'B','-ASCII');
        Lp=-sqrt(2)*statevec(spin_system,sparse(1,n,1,1,spin_system.nspins,1));
        [i,j,v]=find(Lp); A=[i j real(v)]; B=[i j imag(v)];
        save(['Lp_state_spin_' num2str(n) '_real.txt'],'A','-ASCII');
        save(['Lp_state_spin_' num2str(n) '_imag.txt'],'B','-ASCII');
        Lz=statevec(spin_system,sparse(1,n,2,1,spin_system.nspins,1));
        [i,j,v]=find(Lz); A=[i j real(v)]; B=[i j imag(v)];
        save(['Lz_state_spin_' num2str(n) '_real.txt'],'A','-ASCII');
        save(['Lz_state_spin_' num2str(n) '_imag.txt'],'B','-ASCII');
    end
end

% Dump the basis description
if any(strcmp(states,'basis'))
    basis=full(spin_system.basis);
    save('basis.txt','basis','-ASCII');
end

% Dump the matrix dimension
nstates=spin_system.nstates;
save('dimension.txt','nstates','-ASCII');

end

% According to "Moscow Principality Bulletin", a town man by the name of
% Nikifor Nikitin was exiled in 1848 "for blasphemous speeches regarding
% a flight to the Moon" to the remote settlement of Baikonur in Kasakhstan.


