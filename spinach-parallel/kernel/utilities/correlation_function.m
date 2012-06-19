% Rotational correlation function, normalized to be the correlation
% function between second rank Wigner functions.
%
% ilya.kuprov@oerc.ox.ac.uk

function  G=correlation_function(spin_system,k,m,p,q,tau)

% Select the rotational diffusion model
switch numel(spin_system.rlx.tau_c)
    
    case 1
        
        % Use the isotropic rotational diffusion model
        D=1/(6*spin_system.rlx.tau_c);
        G=(1/5)*krondelta(k,p)*krondelta(m,q)*exp(-6*D*tau);
        
    case 2
        
        % Use the axial rotational diffusion model
        D_ax=1/(6*spin_system.rlx.tau_c(1)); D_eq=1/(6*spin_system.rlx.tau_c(2));
        G=(1/5)*krondelta(k,p)*krondelta(m,q)*exp(-(6*D_eq+((k-3)^2)*(D_ax-D_eq))*tau);
        
    case 3
        
        % Use the anisotropic rotational diffusion model
        Dxx=1/(6*spin_system.rlx.tau_c(1));
        Dyy=1/(6*spin_system.rlx.tau_c(2));
        Dzz=1/(6*spin_system.rlx.tau_c(3));
        
        % Refuse to process degenerate cases
        if (abs(Dxx-Dyy)<1e-6*mean([Dxx Dyy Dzz]))||...
           (abs(Dyy-Dzz)<1e-6*mean([Dxx Dyy Dzz]))||...
           (abs(Dzz-Dxx)<1e-6*mean([Dxx Dyy Dzz]))
            error('correlation_function: the three rotational correlation times must be different.');
        end
        
        % Compute the decay rates
        delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz);
        lambda(1)=4*Dxx+Dyy+Dzz;
        lambda(2)=Dxx+4*Dyy+Dzz;
        lambda(3)=Dxx+Dyy+4*Dzz;
        lambda(4)=2*Dxx+2*Dyy+2*Dzz-2*delta;
        lambda(5)=2*Dxx+2*Dyy+2*Dzz+2*delta;
        lambda_p=sqrt(2/3)*(Dxx+Dyy-2*Dzz+2*delta)/(Dxx-Dyy);
        lambda_m=sqrt(2/3)*(Dxx+Dyy-2*Dzz-2*delta)/(Dxx-Dyy);
        
        % Compute the coefficients
        a(1,2)=1/sqrt(2); a(1,4)=1/sqrt(2);
        a(2,2)=-1/sqrt(2); a(2,4)=1/sqrt(2);
        a(3,1)=-1/sqrt(2); a(3,5)=1/sqrt(2);
        a(4,1)=1/sqrt(2+lambda_m^2); a(4,3)=lambda_m/sqrt(2+lambda_m^2); a(4,5)=1/sqrt(2+lambda_m^2);
        a(5,1)=1/sqrt(2+lambda_p^2); a(5,3)=lambda_p/sqrt(2+lambda_p^2); a(5,5)=1/sqrt(2+lambda_p^2);
        
        % Compute the correlation function
        G=0;
        for j=1:5
            G=G+(1/5)*krondelta(m,q)*a(j,k)*a(j,p)*exp(-lambda(j)*tau);
        end

end

% “But they are useless. They can only give you answers.”
%
% Pablo Picasso, about computers.

