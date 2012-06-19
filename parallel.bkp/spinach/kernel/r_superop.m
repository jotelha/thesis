% The relaxation superoperator. Computes the Redfield superoperator by
% numerical evaluation of the integral in the Bloch-Redfield-Wangsness
% master equation in Liouville space. All options are set during the 
% call to create.m function.
%
% The function has been written for minimal memory footprint and performs
% aggressive memory recycling. Relaxation superoperator dimensions in
% excess of 1,000,000 are possible.
%
% ilya.kuprov@oerc.ox.ac.uk
% luke.edwards@chem.ox.ac.uk

function R=r_superop(spin_system)

% Click forward for output
spin_system=click(spin_system,'forward');

% If the correlation time in Redfield theory is set to zero, return a zero matrix
if strcmp(spin_system.rlx.theory,'redfield')
    if (~isfield(spin_system.rlx,'tau_c'))||(norm(spin_system.rlx.tau_c)==0)
        spin_system.rlx.theory='none';
    end
end

% Compute the relaxation superoperator
switch spin_system.rlx.theory

    case 'none'
        
        % Update the user
        report(spin_system,'r_superop: relaxation superoperator set to zero.');
        
        % Just initialize to zero
        R=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,0);
        
    case 'damp'
        
        % Update the user
        report(spin_system,['r_superop: non-selective damping at ' num2str(spin_system.rlx.damp_rate) ' Hz for all states except identity.']);
        
        % Damp everything except the unit state
        R=-spin_system.rlx.damp_rate*speye(spin_system.bas.nstates); R(1,1)=0;
        
    case 't1_t2'
        
        % Update the user
        report(spin_system,'r_superop: extended T1,T2 approximation with user-supplied rates for all states except identity.');
        
        % Prealocate the relaxation rate array
        relaxation_rates=zeros(spin_system.bas.nstates,1);
        
        % Compute the ranks and projections
        [L,M]=lin2lm(spin_system.bas.basis);
        
        % Inspect every state and assign the relaxation rate
        for n=1:spin_system.bas.nstates
            
            % Initialize the relaxation rate to zero
            current_rlx_rate=0;
            
            % Loop over the constituent spins
            for k=1:spin_system.comp.nspins
                % Inspect state rank and projection
                if L(n,k)==0
                    % Spins in identity states do not contribute
                else
                    % Longitudinal states contribute R1, transverse states contribute R2
                    switch M(n,k)
                        case 0
                            current_rlx_rate=current_rlx_rate+spin_system.rlx.r1(k);
                        otherwise
                            current_rlx_rate=current_rlx_rate+spin_system.rlx.r2(k);
                    end
                end
            end
            relaxation_rates(n)=current_rlx_rate;
            
        end
        
        % Deallocate the rank and projection arrays
        clear('L','M');
        
        % Form the relaxation superoperator
        R=spdiags(-relaxation_rates,0,spin_system.bas.nstates,spin_system.bas.nstates);

    case 'redfield'
        
        % Preallocate the relaxation superoperator
        R=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,10*spin_system.bas.nstates);
        
        % Get the rotational basis, including the non-secular terms
        report(spin_system,'r_superop: computing the lab frame Hamiltonian superoperator...');
        [L0,Q]=h_superop(secularity(spin_system,'keep_all'));
        
        % Determine the largest and the smallest correlation time
        tau_c_max=max(spin_system.rlx.tau_c);
        tau_c_min=min(spin_system.rlx.tau_c);
        
        % Determine the norm of L0
        L0_norm=norm(L0,'inf');
        
        % Set the upper integration limit according to the accuracy goal
        upper_limit=tau_c_max*log(1/spin_system.tols.rlx_integration);
        
        % Set the number of integration steps according to the accuracy goal
        nsteps=ceil((1/(10*(spin_system.tols.rlx_integration^(1/6))))*...
                   ((tau_c_max*log(1/spin_system.tols.rlx_integration)/min([tau_c_max 1/L0_norm]))^(7/6)));
        
        % Get the numerical integration step
        timestep=upper_limit/nsteps;
        
        % Kill the terms in the static Liouvillian that are irrelevant on the time scale of the integration
        L0=L0.*(upper_limit*abs(L0)>1);
        
        % Report back to the user
        report(spin_system,['r_superop: infinity-norm of the static Hamiltonian superoperator ' num2str(L0_norm)]);
        report(spin_system,['r_superop: numerical integration time step ' num2str(timestep) ' seconds.']);
        report(spin_system,['r_superop: ' num2str(nsteps) ' integration steps will be taken.']);
        
        % Get the static Liouvillian propagator
        P=propagator(spin_system,L0,timestep/4);
        
        % Deallocate the static Liouvillian
        clear('L0');
        
        % Loop over the rotational basis operators
        for k=1:5
            for m=1:5
                for p=1:5
                    for q=1:5
                        
                        % Pick the correlation function
                        G=@(tau)correlation_function(spin_system,k,m,p,q,tau);
                        
                        % Compute the term in the rotational expansion sum
                        if significant(Q{k,m},'operator',eps)&&...
                           significant(Q{p,q},'operator',eps)&&...
                           significant(G(tau_c_min),'scalar',eps)
                    
                            % Preallocate the integration result
                            A=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,10*spin_system.bas.nstates);
                    
                            % Set the initial value for the the operator between
                            % the exponentials in the BRW integral
                            B=Q{p,q}';
                            
                            % Compute the BRW integral using O(h^7) Boole's quadrature
                            for n=0:(nsteps-1)
                                A=A+((timestep/90)* 7*G((n+0.00)*timestep))*B;
                                B=clean_up(spin_system,P'*B*P,spin_system.tols.liouv_zero);
                                A=A+((timestep/90)*32*G((n+0.25)*timestep))*B;
                                B=clean_up(spin_system,P'*B*P,spin_system.tols.liouv_zero);
                                A=A+((timestep/90)*12*G((n+0.50)*timestep))*B;
                                B=clean_up(spin_system,P'*B*P,spin_system.tols.liouv_zero);
                                A=A+((timestep/90)*32*G((n+0.75)*timestep))*B;
                                B=clean_up(spin_system,P'*B*P,spin_system.tols.liouv_zero);
                                A=A+((timestep/90)* 7*G((n+1.00)*timestep))*B;
                            end
                    
                            % Deallocate unnecessary variables
                            clear('B');
                            
                            % Add the result to the total
                            R=R-clean_up(spin_system,Q{k,m}*A,spin_system.tols.liouv_zero);
                    
                        end
                        
                    end
                end
            end
        end
        
        % Deallocate variables
        clear('Q','A','B','P');
        
    otherwise
        error('r_superop: unknown relaxation superoperator type.');
        
end

% Print matrix density statistics
if ~strcmp(spin_system.rlx.theory,'none')
    report(spin_system,['r_superop: full relaxation superoperator density ' num2str(100*nnz(R)/numel(R)) '%, nnz(R)=' num2str(nnz(R))]);
end

% Decide the fate of the dynamic frequency shifts
if strcmp(spin_system.rlx.theory,'redfield')
    switch spin_system.rlx.dfs
        case 'keep'
            report(spin_system,'r_superop: dynamic frequency shifts have been kept.');
        case 'ignore'
            R=real(R);
            report(spin_system,'r_superop: dynamic frequency shifts have been ignored.');
        otherwise
            error('r_superop: invalid value of the inter.rlx_dfs parameter.');
    end
end

% Decide the fate of the off-diagonal components
if strcmp(spin_system.rlx.theory,'redfield')
    
    switch spin_system.rlx.keep
        case 'diagonal'
            
            % Pull out the diagonal
            R=diag(diag(R));
            
            % Inform the user
            report(spin_system,'r_superop: all cross-relaxation terms have been ignored.');
            
        case 'kite'
            
            % Pull out the diagonal
            R_diag=diag(diag(R)); R_offdiag=R-R_diag;
            
            % Compile the index of all longitudinal spin orders
            [~,M]=lin2lm(spin_system.bas.basis);
            longitudinals=double(sparse(sum(abs(M),2)==0));
            
            % Apply the mask in a memory-friendly way
            for n=1:size(R,2)
                R_offdiag(:,n)=R_offdiag(:,n).*(longitudinals*longitudinals(n));
            end
            
            % Reassemble the superoperator
            R=R_diag+R_offdiag;
            
            % Inform the user
            report(spin_system,'r_superop: transverse cross-relaxation terms have been ignored.');
            
        case 'secular'
            
            % Compute state carrier frequencies
            [~,M]=lin2lm(spin_system.bas.basis);
            frequencies=sum(repmat(spin_system.inter.basefrq,spin_system.bas.nstates,1).*M,2);
            
            % Apply the secularity mask in a memory-friendly way
            for omega=unique(frequencies)'
                mask=sparse(abs(frequencies-omega)<1e-6);
                for n=find(mask)'
                    R(:,n)=R(:,n).*mask; %#ok<SPRIX>
                end
            end
            
            % Inform the user
            report(spin_system,'r_superop: non-secular cross-relaxation terms have been ignored.');
            
        case 'full'
            
            % Inform the user
            report(spin_system,'r_superop: returning full relaxation superoperator (lab frame simulations only).');
            
    end
    
    % Print matrix sparsity statistics
    report(spin_system,['r_superop: final relaxation superoperator density ' num2str(100*nnz(R)/numel(R)) '%, nnz(R)=' num2str(nnz(R))]);   
    
end

% Decide the equilibrium state
if ~strcmp(spin_system.rlx.theory,'none')
    switch spin_system.rlx.equilibrium
        case 'zero'
            report(spin_system,'r_superop: WARNING -- the spin system will relax to the all-zero state.');
        case 'thermal'
            report(spin_system,'r_superop: the system will relax to the thermal equilibrium.');
            R(:,1)=-R*equilibrium(spin_system);
        otherwise
            error('r_superop: unknown equilibrium specification.');
    end
end

% Set the absorptive boundary if required (undocumented)
if isfield(spin_system.rlx,'boundary')&&strcmp(spin_system.rlx.boundary,'absorptive')
    
    % Determine the highest spin order in the system
    spin_orders=sum(logical(spin_system.bas.basis),2);
    highest_spin_order=max(spin_orders);
    
    % Find all such spin orders
    boundary_orders=(spin_orders==highest_spin_order);
    
    % Add damping terms to the relaxation superoperator
    R=R-spin_system.rlx.boundary_damping_rate*spdiags(boundary_orders,0,length(boundary_orders),length(boundary_orders));
    
    % Report back to the user
    report(spin_system,['r_superop: ' num2str(spin_system.rlx.boundary_damping_rate) ' Hz damping applied to all ' num2str(highest_spin_order) '-spin orders.']);
    
end

end

% There are horrible people who, instead of solving a problem, tangle it up
% and make it harder to solve for anyone who wants to deal with it. Whoever
% does not know how to hit the nail on the head should be asked not to hit
% it at all.
%
% Friedrich Nietzsche


