% Calculates Clebsch-Gordan coefficients.
%                                
% If physically inadmissible indices are supplied,
% a zero is returned.
%
% A very considerable amount of thought has been given
% to the accuracy and performance of this function.
%
% ilya.kuprov@oerc.ox.ac.uk

function cg=clebsch_gordan(L,M,L1,M1,L2,M2)

% Check the input
grumble(L,M,L1,M1,L2,M2);

% Run a fast bypass for L1=2 (used by the SLE module)
if L1==2
    cg=cg_bypass(L,M,L1,M1,L2,M2); return
end

% Set the default answer
cg=0;

% Match the notation to Varshalovich, Section 8.2.1
c=L; gam=M; a=L1; alp=M1; b=L2; bet=M2;

% Run zero tests (Stage I)
prefactor_is_nonzero=(a+alp>=0)&&(a-alp>=0)&&...
                     (b+bet>=0)&&(b-bet>=0)&&...
                     (c+gam>=0)&&(c-gam>=0)&&...
                      krondelta(gam,alp+bet);
    
% Proceed if appropriate
if prefactor_is_nonzero
    
    % Run zero tests (Stage II)
    delta_is_nonzero=(a+b-c>=0)&&(a-b+c>=0)&&...
                     (-a+b+c>=0)&&(a+b+c+1>=0);
                     
    % Proceed if appropriate
    if delta_is_nonzero
            
        % Run zero tests (Stage III)
        lower_sum_limit=max([alp-a, b+gam-a, 0]);
        upper_sum_limit=min([c+b+alp, c+b-a, c+gam]);
        sum_is_nonzero=(upper_sum_limit>=lower_sum_limit);
    
        % Proceed if appropriate
        if sum_is_nonzero
            
            % Compile a look-up table of gamma functions
            gamma(a+b+c+2)=java.math.BigInteger.ONE; gamma(:)=gamma(a+b+c+2);
            for k=3:(a+b+c+2)
                gamma(k)=gamma(k-1).multiply(java.math.BigInteger.valueOf(k-1));
            end
            
            % Compute prefactor numerator
            numer=java.math.BigInteger.valueOf(2*c+1);
            numer=numer.multiply(gamma(a-b+c+1));
            numer=numer.multiply(gamma(-a+b+c+1));
            numer=numer.multiply(gamma(c+gam+1));
            numer=numer.multiply(gamma(c-gam+1));
            numer=numer.multiply(gamma(a+b-c+1));
            
            % Compute prefactor denominator
            denom=gamma(a+b+c+2);
            denom=denom.multiply(gamma(a+alp+1));
            denom=denom.multiply(gamma(a-alp+1));
            denom=denom.multiply(gamma(b+bet+1));
            denom=denom.multiply(gamma(b-bet+1));
            
            % Perform floating-point division
            accur=length(denom.toString)-length(numer.toString)+64;
            numer=java.math.BigDecimal(numer); denom=java.math.BigDecimal(denom);
            prefactor=numer.divide(denom,accur,java.math.RoundingMode.HALF_UP);
            
            % Compute the sum
            z_sum=java.math.BigDecimal.ZERO;
            for z=lower_sum_limit:upper_sum_limit
                
                % Compute term numerator
                numer=java.math.BigInteger.valueOf((-1)^(b+bet+z));
                numer=numer.multiply(gamma(c+b+alp-z+1));
                numer=numer.multiply(gamma(a-alp+z+1));
                
                % Compute term denominator
                denom=gamma(z+1);
                denom=denom.multiply(gamma(c-a+b-z+1));
                denom=denom.multiply(gamma(c+gam-z+1));
                denom=denom.multiply(gamma(a-b-gam+z+1));
                
                % Perform floating-point division
                numer=java.math.BigDecimal(numer); denom=java.math.BigDecimal(denom);
                sum_term=numer.divide(denom,accur,java.math.RoundingMode.HALF_UP);
                
                % Add the term to the total
                z_sum=z_sum.add(sum_term);
                
            end
            
            % Return to double precision
            cg_sq=z_sum.multiply(z_sum).multiply(prefactor);
            cg=z_sum.signum*sqrt(double(cg_sq));
            
        end
            
    end
        
end

end

% Input checker
function grumble(L,M,L1,M1,L2,M2)

if (double(int64(2*L+1))~=(2*L+1))||(double(int64(2*M+1))~=(2*M+1))||...
   (double(int64(2*L1+1))~=(2*L1+1))||(double(int64(2*M1+1))~=(2*M1+1))||...
   (double(int64(2*L2+1))~=(2*L2+1))||(double(int64(2*M2+1))~=(2*M2+1))
    error('clebsch_gordan: the arguments must be integer or half-integer.');
end

end

% Explicit formulae for the L1=2 case
function cg=cg_bypass(L,M,~,M1,L2,M2)

% Enumerate all cases explicitly
switch M1
    case -2
        if M~=M2-2, cg=0; else
            switch L
                case L2-2
                    cg=sqrt(((-3+L2+M2)*(-2+L2+M2)*(-1+L2+M2)*(L2+M2))/...
                           (L2*(1-L2-4*power(L2,2)+4*power(L2,3))))/2;
                case L2-1
                    cg=-(sqrt(((1+L2-M2)*(-2+L2+M2)*(-1+L2+M2)*(L2+M2))/....
                        (L2*(-1-2*L2+power(L2,2)+2*power(L2,3))))/sqrt(2));
                case L2
                    cg=sqrt(1.5)*sqrt(((1+L2-M2)*(2+L2-M2)*(-1+L2+M2)*(L2+M2))/...
                                     (L2*(-3+L2*(1+4*L2*(2+L2)))));
                case L2+1
                    cg=-(sqrt(((1+L2-M2)*(2+L2-M2)*(3+L2-M2)*(L2+M2))/...
                             (L2*(2+7*L2+7*power(L2,2)+2*power(L2,3))))/sqrt(2));
                case L2+2
                    cg=(power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(2+L2-M2)*(3+L2-M2)*(4+L2-M2))/...
                             (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4))))/2;
                otherwise
                    cg=0;
            end
        end
    case -1
        if M~=M2-1, cg=0; else
            switch L
                case L2-2
                    cg=-sqrt(((L2-M2)*(-2+L2+M2)*(-1+L2+M2)*(L2+M2))/...
                            (L2*(1-L2-4*power(L2,2)+4*power(L2,3))));
                case L2-1
                    cg=((1+L2-2*M2)*sqrt(((-1+L2+M2)*(L2+M2))/...
                                        (L2*(-1-2*L2+power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2
                    cg=sqrt(1.5)*sqrt(((1+L2-M2)*(L2+M2))/...
                                     (L2*(-3+L2*(1+4*L2*(2+L2)))))*(-1+2*M2);
                case L2+1
                    cg=-((power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(2+L2-M2))/...
                               (L2*(1+L2)*(2+L2)*(1+2*L2)))*(L2+2*M2))/sqrt(2));
                case L2+2
                    cg=power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(2+L2-M2)*(3+L2-M2)*(1+L2+M2))/...
                            (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4)));
                otherwise
                    cg=0;
            end
        end
    case 0
        if M~=M2, cg=0; else
            switch L
                case L2-2
                    cg=sqrt(1.5)*sqrt(((-1+L2-M2)*(L2-M2)*(-1+L2+M2)*(L2+M2))/...
                                     (L2*(1-L2-4*power(L2,2)+4*power(L2,3))));
                case L2-1
                    cg=sqrt(3)*M2*sqrt((power(L2,2)-power(M2,2))/...
                                      (L2*(-1-2*L2+power(L2,2)+2*power(L2,3))));
                case L2
                    cg=(-(L2*(1+L2))+3*power(M2,2))/...
                                       sqrt(L2*(-3+L2*(1+4*L2*(2+L2))));
                case L2+1
                    cg=-(power(-1,-2*L2+2*M2)*sqrt(3)*M2*sqrt(((1+L2-M2)*(1+L2+M2))/...
                                                             (L2*(1+L2)*(2+L2)*(1+2*L2))));
                case L2+2
                    cg=power(-1,-2*L2+2*M2)*sqrt(3)*sqrt(((1+L2-M2)*(2+L2-M2)*(1+L2+M2)*(2+L2+M2))/...
                            (12+50*L2+70*power(L2,2)+40*power(L2,3)+8*power(L2,4)));
                otherwise
                    cg=0;
            end
        end
    case 1
        if M~=M2+1, cg=0; else
            switch L
                case L2-2
                    cg=-sqrt(((-2+L2-M2)*(-1+L2-M2)*(L2-M2)*(L2+M2))/...
                            (L2*(1-L2-4*power(L2,2)+4*power(L2,3))));
                case L2-1
                    cg=-((sqrt(((-1+L2-M2)*(L2-M2))/...
                        ((-1+L2)*L2*(1+L2)*(1+2*L2)))*(1+L2+2*M2))/sqrt(2));
                case L2
                    cg=-(power(-1,-2*L2+2*M2)*sqrt(1.5)*sqrt(((L2-M2)*(1+L2+M2))/...
                              (L2*(-3+L2*(1+4*L2*(2+L2)))))*(1+2*M2));
                case L2+1
                    cg=(power(-1,-2*L2+2*M2)*(L2-2*M2)*sqrt(((1+L2+M2)*(2+L2+M2))/...
                             (L2*(2+7*L2+7*power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2+2
                    cg=power(-1,-2*L2+2*M2)*sqrt(((1+L2-M2)*(1+L2+M2)*(2+L2+M2)*(3+L2+M2))/...
                            (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4)));
                otherwise
                    cg=0;
            end
        end
    case 2
        if M~=M2+2, cg=0; else
            switch L
                case L2-2
                    cg=sqrt(((-3+L2-M2)*(-2+L2-M2)*(-1+L2-M2)*(L2-M2))/...
                           (L2*(1-L2-4*power(L2,2)+4*power(L2,3))))/2;
                case L2-1
                    cg=(power(-1,-2*L2+2*M2)*sqrt(((-2+L2-M2)*(-1+L2-M2)*(L2-M2)*(1+L2+M2))/...
                             (L2*(-1-2*L2+power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2
                    cg=power(-1,-2*L2+2*M2)*sqrt(1.5)*sqrt(((-1+L2-M2)*(L2-M2)*(1+L2+M2)*(2+L2+M2))/...
                                                          (L2*(-3+L2*(1+4*L2*(2+L2)))));
                case L2+1
                    cg=(power(-1,-2*L2+2*M2)*sqrt(((L2-M2)*(1+L2+M2)*(2+L2+M2)*(3+L2+M2))/...
                             (L2*(2+7*L2+7*power(L2,2)+2*power(L2,3)))))/sqrt(2);
                case L2+2
                    cg=(power(-1,-2*L2+2*M2)*sqrt(((1+L2+M2)*(2+L2+M2)*(3+L2+M2)*(4+L2+M2))/...
                             (6+25*L2+35*power(L2,2)+20*power(L2,3)+4*power(L2,4))))/2;
                otherwise
                    cg=0;
            end
        end
    otherwise
        cg=0;
end

end

% The miracle of the appropriateness of the language of mathematics for the
% formulation of the laws of physics is a wonderful gift which we neither
% understand nor deserve. We should be grateful for it and hope that it will
% remain valid in future research and that it will extend, for better or for
% worse, to our pleasure, even though perhaps also to our bafflement, to
% wide branches of learning.
%
% Eugene Wigner, http://dx.doi.org/10.1002/cpa.3160130102

