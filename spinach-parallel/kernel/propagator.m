% Calculates exponential propagators and propagator derivatives. A
% three-parameter call:
%
%   P=propagator(spin_system,L,timestep)
%
% returns exp(-i*L*t). A four-parameter call also returns the propagator
% derivatives with respect to Liouvillian parameters:
%
%   [P,dP]=prop_gen(spin_system,L,timestep,dL)
%
% where dL is the cell array of Liouvillian derivatives with respect
% to the parameters in question and dP is a cell array of the corresponding
% derivative propagators.
%
% ilya.kuprov@oerc.ox.ac.uk

function [P,dP]=propagator(spin_system,L,timestep,dL)

% Click forward for output
spin_system=click(spin_system,'forward');

% Set the shorthand for -i*L*dt
A=-1i*L*timestep;

% Set the shorthands for the derivatives
if nargin==4
    B=cell(size(dL)); 
    for n=1:length(dL)
        B{n}=-1i*dL{n}*timestep;
    end
end

% Decide between the small-matrix and the large-matrix paths
if (nargin==3)&&(size(A,1)<spin_system.tols.small_matrix)
    
    % Use Matlab's expm function for small matrices
    P=expm(full(A));

else
    
    % Inform the user about matrix densities
    report(spin_system,['propagator: Liouvillian density: ' num2str(100*nnz(L)/numel(L)) ' %']);
    if nargin==4
        report(spin_system,['propagator: Liouvillian derivative(s) density: ' num2str(100*cellfun(@nnz,dL)./cellfun(@numel,dL)) ' %']);
    end
    
    % Estimate the biggest eigenvalue
    maxeig=normest(A,spin_system.tols.normest_tol);
    
    % Warn the user if the biggest eigenvalue is too big
    if maxeig>1024
        
        % If the user is really pusing it, take precautionary measures
        report(spin_system,'propagator: WARNING - the time step requested greatly exceeds the timescale of system dynamics.');
        report(spin_system,'propagator: WARNING - all exponentiation tolerances will be set to 1e-14.');
        spin_system.tols.prop_norm=1e-14;
        spin_system.tols.prop_chop=1e-14;
        spin_system.tols.prop_zero=1e-14;
        
    elseif maxeig>16
        
        % Inform the user just in case
        report(spin_system,'propagator: WARNING - the time step requested exceeds the timescale of system dynamics.');
        
    end
    
    % Apply scaling if necessary
    if maxeig>1
        
        % Determine the scaling factor
        scaling_factor=2^ceil(log2(maxeig));
        
        % Scale the -iLdt
        A=A/scaling_factor;
        
        % Scale the derivatives
        if nargin==4
            for n=1:length(B), B{n}=B{n}/scaling_factor; end
        end
        
        % Update the user
        report(spin_system,['propagator: Ldt eigenvalue estimate outside monotonic convergence radius. Scaling and squaring applied, scaling_factor=' num2str(scaling_factor)]);
        
    end
    
    % Get the propagator
    switch spin_system.tols.exponentiation
        
        case 'taylor'
    
            % Run the Taylor series procedure for the propagator
            P=speye(size(A)); next_term=speye(size(A)); next_term_norm=1; n=1;
            while next_term_norm > spin_system.tols.prop_norm
                
                % Compute the next term
                next_term=next_term*A/n;
                
                % Run the clean-up procedure
                next_term=clean_up(spin_system,next_term,spin_system.tols.prop_chop);
                
                % Compute the residual norm
                next_term_norm=max(max(abs(next_term)));
                
                % Add to the total and increment the counter
                P=P+next_term; n=n+1;
                
            end
            report(spin_system,['propagator: Taylor series converged in ' num2str(n) ' iterations.']);
    
        case 'chebyshev'
            
            % Start the Chebyshev series
            Tnm1=speye(size(A)); Tn=-1i*A;
            
            % Write the first two terms into the propagator
            P=besselj(0,1)*Tnm1+2i*besselj(1,1)*Tn;
            
            % Initialize the convergence flag and iteration counter
            next_term_norm=1; n=2;
            
            % Get the rest of the terms
            while next_term_norm > spin_system.tols.prop_norm
                
                % Update the Chebyshev polynomials
                Tnp1=-2i*A*Tn-Tnm1; Tnm1=Tn; Tn=Tnp1;
                
                % Run the clean-up procedure
                Tn=clean_up(spin_system,Tn,spin_system.tols.prop_chop);
                
                % Get the next term
                next_term=2*(1i^n)*besselj(n,1)*Tn;
                
                % Compute the residual norm
                next_term_norm=max(max(abs(next_term)));
                
                % Add to the total and increment the counter
                P=P+next_term; n=n+1;
                
            end
            report(spin_system,['propagator: Chebyshev series converged in ' num2str(n) ' iterations.']);
            
        otherwise
            
            error('propagator: unknown matrix exponentiation method.');
            
    end
    
    % Run the Hausdorff series procedure for the derivative propagator(s)
    if nargin==4
        dP=cell(size(B));
        for r=1:length(B)
            
            % Preallocate the array
            dP{r}=spalloc(spin_system.bas.nstates,spin_system.bas.nstates,spin_system.bas.nstates);
            
            % Converge the Hausdorff series
            next_term=B{r}; next_term_norm=1; n=1;
            while next_term_norm > spin_system.tols.derivative_prop_norm
                
                % Add to the total
                dP{r}=dP{r}+next_term; n=n+1;
                
                % Compute the next term
                next_term=(next_term*A-A*next_term)/n;
                
                % Run the clean-up procedure and check the density
                next_term=clean_up(spin_system,next_term,spin_system.tols.derivative_prop_chop);
                
                % Compute the residual norm
                next_term_norm=max(abs(next_term(:)));
                
            end
            
            % Compute the derivatives
            dP{r}=P*dP{r};
            
            % Run the clean-up procedure
            dP{r}=clean_up(spin_system,dP{r},spin_system.tols.derivative_prop_chop);
            
            % Inform the user
            report(spin_system,['propagator: derivative propagator # ' num2str(r) ' converged in ' num2str(n) ' iterations.']);
            
        end
    end
    
    % Square back if necessary
    if maxeig>1
        
        % Update the user
        report(spin_system,'propagator: squaring the propagator up to the original time step...');
        
        % Process the propagators
        for n=1:ceil(log2(maxeig))
            
            % Process the derivative propagators
            if nargin==4
                for k=1:length(dP)
                    
                    % Run the squaring step
                    dP{k}=dP{k}*P+P*dP{k};
                    
                    % Run the clean-up procedure
                    dP{r}=clean_up(spin_system,dP{r},spin_system.tols.derivative_prop_chop);
                    
                end
            end
            
            % Run the squaring step
            P=P^2;
            
            % Run the clean-up procedure
            P=clean_up(spin_system,P,spin_system.tols.prop_chop);
            
        end
        
    end
    
    % Do the final clean-up and reporting on the main propagator
    P=clean_up(spin_system,P,spin_system.tols.prop_zero);
    report(spin_system,['propagator: propagator density: ' num2str(100*nnz(P)/numel(P)) ' %']);
    
    % Do the final clean-up and reporting on the derivative propagator(s)
    if nargin==4
        for n=1:length(dP)
            dP{r}=clean_up(spin_system,dP{r},spin_system.tols.derivative_prop_zero);
            report(spin_system,['propagator: derivative propagator # ' num2str(n) ' density: ' num2str(100*nnz(dP{r})/numel(dP{r})) ' %']);
        end
    end
    
end

end

% To preserve one's mind intact through a modern college education is a
% test of courage and endurance, but the battle is worth it and the stakes
% are the highest possible to man: the survival of reason.
%
% Ayn Rand

