% Implementation of the Direct Expectation values via Chebyshev (DEC)
% algorithm for propagation based upon code by G. Mazzi. NB. unlikely to
% converge for non-Hermitian L.
%
% Inputs:   spin_system - /Spinach/ struct
%           L - Liouvillian
%           coil - detection state
%           rho - initial state
%           timestep - timetep required
%           nsteps - number of steps required
%
% luke.edwards@chem.ox.ac.uk
% giacomo.mazzi@cs.kuleuven.be

function fid=dec(spin_system,L,coil,rho,timestep,nsteps)

% Click forward for output
spin_system=click(spin_system,'forward');

% Time axis
T=0:timestep:nsteps*timestep;

% Preallocate return array
fid=zeros(1,length(T));

% Scale Liouvillian
report(spin_system,'dec: scaling the Liouvillian...')
[S,D,Ls]=scale_L(L,0);

% First point is ready
fid(1)=coil'*rho;

% Generate R, and the expectation value for the final time.
report(spin_system,'dec: calculating Chebyshev polynomial for final time...')
[fid(end),R,n]=chebfinal(Ls,coil,rho,T(end)*D,eps,0);
report(spin_system,['dec: requires ', num2str(n-1), ' Chebyshev steps.']);

% Sanity check on output
if norm(R(1:floor(end/2)))<10*norm(R(ceil(end/2):end))
    report(spin_system,'dec: WARNING: check eigenvalues contained in ellipse.')
end

% Reevaluate ck for all previous times and calculate expectation values.
report(spin_system,'dec: interpolating earlier times...')
fid(2:end-1)=chebterms(T*D,R,eps);

% Multiply by exponential scaling factor
fid=fid.*exp(-1i*S*T);

end

%%
function [S,D,Ls]=scale_L(L,debug_flag)
% Scale the Liouvillian, L, such that the series converges. Non-Hermitian L
% are scaled using a more general formalism, c.f. Huang, et al.

eigenvalue=eigs(L,1,'lm');

if abs(imag(eigenvalue))>eps 
    try
        % Because the real part of the eigenvalues are often so much bigger
        % than the imaginary, this often fails to converge
        S=0.7i*imag(eigs(L,1,'si'));
        real_eig=1.2*real(eigenvalue);

    catch exception
        try
            % Rescale matrix and try again
            L_im=real(L)./real(eigenvalue)+1i*imag(L);
            S=0.5i*imag(eigs(L_im,1,'si'));
            real_eig=1.3*real(eigenvalue);
            disp(['WARNING: ' exception.identifier])
            disp(exception.message)
            
        catch exception2
            
            % Be very conservative with eigenvalue estimates 
            S=10i*imag(eigenvalue);
            real_eig=2*real(eigenvalue);
            disp(['WARNING: ' exception2.identifier])
            disp(exception.message)
        end
    end
    
    % Half the distance between the two foci of the ellipse
    D=sqrt(real_eig^2-imag(S)^2);
    
else
    
    % Matrix is Hermitian
    real_eig=1.01*eigenvalue;
    S=0;
    D=real_eig;
    
end

% Plot the bounding ellipse and the maximal eigenvalues for debugging
% purposes
if exist('debug_flag','var')&&debug_flag
    figure, plot(eigs(L,min([250,length(L)-1])),'x'), hold on
    t=linspace(0,2*pi,100);
    x=(real_eig)*cos(t); y=(imag(S))*sin(t);
    plot(x,y+imag(S))
    hold off
    drawnow, disp('Press any key to continue...'), pause
end

% Scale matrix
Ls=(L-S*speye(size(L)))/D;

end

%%
function [fid,R,n]=chebfinal(Ls,coil,rho,t,tol,debug_flag)
% Generates the array of complex doubles that are used to generate the fid
% for all subsequent times.
%
% Inputs:   tol - tolerance for the convergence
%           Ls - scaled Liouvillian
%           coil - detection state
%           rho - initial density matrix
%           t - final time
%           R - array containing <coil|Tn|rho>
%
% Outputs:  fid - final point of the fid
%           R - see above
%           n - number of steps taken for convergence

if ~exist('tol','var'), tol=eps; end

% Evaluate T_0 and T_1
Tnm1=rho;  %0 term for cheby
Tn=Ls*rho; %1 term for cheby

% Generate initial Bessel functions
Jnm1=besselj(0,t);
Jn=besselj(1,t);

% Initialize R (N.B. R likely to be larger than this)
R=zeros(1,ceil(abs(t)));

% First two R_k are ready
R(1)=coil'*Tnm1;
R(2)=-2i*(coil'*Tn);

% Update fid with first two terms
fid=Jnm1*R(1)+Jn*R(2);

% Convergence of Bessel functions exponential when n > t, so only check for
% convergence after that. Avoid using recurrance relation for Bessel
% functions after n>floor(abs(t)) because it becomes unstable.
n=1;
t_rec=1/t;
for n=2:floor(abs(t))
    
    % Update for recurrence
    Jnm2=Jnm1;
    Jnm1=Jn;
    
    % Recurrence relation for Bessel functions
    Jn=2*(n-1)*t_rec*Jnm1-Jnm2;
    
    % Recurrence relation for Chebyshev polynomials
    Tnp1=2*(Ls*Tn)-Tnm1;
    
    % Store expectation value
    R(n+1)=2*(-1i)^n*(coil'*Tnp1);
    
    % Update Tn and Tnm1
    Tnm1=Tn;
    Tn=Tnp1;
    
    % Final fid point updated
    fid=fid+Jn*R(n+1);
    
end

ck2=norm(Jn)+norm(Jnm1);
while (n<ceil(abs(t)))||(ck2>tol)
    
    n=n+1;
    
    % Update for norm-check
    Jnm1=Jn;
    
    % Non-recursive as unstable
    Jn=besselj(n,t);
    
    % Update norm-check variable
    ck2=norm(Jn)+norm(Jnm1);
    
    % Recurrence for Chebyshev polynomials
    Tnp1=2*(Ls*Tn)-Tnm1;
    
    % Store expectation value
    R(n+1)=2*(-1i)^n*(coil'*Tnp1);
    
    % Update Tn and Tnm1
    Tnm1=Tn;
    Tn=Tnp1;
    
    % Update final fid point
    fid=fid+Jn*R(n+1);
    
end

if exist('debug_flag','var')&&debug_flag
    plot([real(R);imag(R)].'), drawnow
    disp('Press any key to continue...'), pause
end

end


%%
function [fid]=chebterms(T,R,tol)
% Utilises the vector, R, generated by chebfinal() to form the fid at a
% time, t < t_final. Runs evaluation backwards to avoid numerical
% instabilities.

if ~exist('tol','var'), tol=1e-6; end
fid=zeros(1,length(T)-2);

% Almost certainly not particularly well load-balanced.
parfor j=1:length(T)-2
    
    t=T(j+1);
    
    % Estimate number of Chebyshev terms required
    N=min([ceil((1.01)*abs(t))+40,length(R)-1]);
    
    % If guessed term too small then can't begin recurrence; calculate
    % lower Bessel functions until one is large enough. Go in steps of 2 
    % to increase speed
    JN=besselj(N-2,t);
    Jprev=besselj(N,t);
    N=N-2;
    while abs(JN)<tol
        N=N-2;
        Jprev=JN;
        JN=besselj(N,t);
    end
    
    % Interpolate the Bessel function skipped by going in steps of 2
    J=zeros(1,N+3);
    J(end)=Jprev;
    J(end-1)=0.5*t/(N+1)*(JN+Jprev);
    J(end-2)=JN;
    
    t_rec=1/t;
    for n=N:-1:1
        J(n)=2*n*t_rec*J(n+1)-J(n+2);
    end
    
    % take vectorized product; multiplication in for loop too slow
    fid(j)=fid(j)+sum(J.*R(1:N+3));
    
end

end
