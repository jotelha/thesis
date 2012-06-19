% Performs fid apodization.
%
% ilya.kuprov@oerc.ox.ac.uk

function fid=apodization(fid,window_type,parameters)

switch window_type
    
    case 'none-1d'
        fid(1)=fid(1)/2;
    
    case 'crisp-1d'
        fid(1)=fid(1)/2;
        fid=fid.*cos(linspace(0,pi/2,length(fid))').^8;
    
    case 'exponential-1d'
        fid(1)=fid(1)/2;
        fid=fid.*exp(-parameters(1)*linspace(0,1,length(fid))');
    
    case 'gaussian-1d'
        fid(1)=fid(1)/2;
        fid=fid.*exp(-parameters(1)*(linspace(0,1,length(fid))').^2);
        
    case 'sinebell-1d'
        fid(1)=fid(1)/2;
        fid=fid.*cos(linspace(0,pi/2,length(fid))');
        
    case 'kaiser-1d'
        fid(1)=fid(1)/2;
        fid=fid.*kaiser(length(fid),parameters(1));
        
    case 'hamming-1d'
        fid(1)=fid(1)/2;
        fid=fid.*hamming(length(fid));
        
    case 'none-2d'
        fid(1,1)=fid(1,1)/2;
        
    case 'crisp-2d'
        fid(1,1)=fid(1,1)/2;
        [ramp_f1,ramp_f2]=meshgrid(linspace(1,0,size(fid,2)),linspace(1,0,size(fid,1)));
        fid=fid.*sin((pi/2)*ramp_f1.*ramp_f2).^8;
    
    case 'exponential-2d'
        fid(1,1)=fid(1,1)/2;
        decay1=exp(-parameters(1)*linspace(0,1,size(fid,2)));
        decay2=exp(-parameters(1)*linspace(0,1,size(fid,1)));
        fid=fid.*kron(decay1,decay2');
    
    case 'gaussian-2d'
        fid(1,1)=fid(1,1)/2;
        decay1=exp(-parameters(1)*linspace(0,1,size(fid,2)).^2);
        decay2=exp(-parameters(1)*linspace(0,1,size(fid,1)).^2);
        fid=fid.*kron(decay1,decay2');
        
    case 'sinebell-2d'
        fid(1,1)=fid(1,1)/2;
        [ramp_f1,ramp_f2]=meshgrid(linspace(1,0,size(fid,2)),linspace(1,0,size(fid,1)));
        fid=fid.*sin((pi/2)*ramp_f1.*ramp_f2);
        
    otherwise
        error(['apodization: function ' window_type ' not implemented.']);

end

end

% I have had my results for a long time: but I do not yet know how I am to arrive at them.
% 
% Carl Friedrich Gauss

