% Reads wavegauss.pk files. Arguments:
%
%              filename   -   a string containing the name of the file
%              npoints    -   if specified, harmonic downsampling would 
%                             be performed to this number of points.
% 
% The 'integral' output refers to the integral of the RF amplitude across
% the pulse. This may be useful to estimate your calibration, but for
% complicated shaped pulses it must be done 'the spectrometer way'.
%
% Konstantin Pervushin (p@trosy.com)
% Russell Bertrand     (s070062@ntu.edu.sg)
% Ilya Kuprov          (ilya.kuprov@durham.ac.uk)

function [amplitudes,phases,integral]=read_wave(filename,npoints)

% Correct the slashes to match the system architecture
if ispc
    filename=regexprep(filename,'/','\');
else
    filename=regexprep(filename,'\','/');
end

% Open the file
wavefile = fopen(filename, 'r');
if wavefile<0
    error(['read_wave: ' filename ' does not exist or could not be opened.']);
else
    disp(['read_wave: opened waveform ' filename]);
end;

% Read the file
textscan(wavefile, '#%s','commentStyle','#');
waveform = textscan(wavefile, '%f, %f');
fclose(wavefile);

% Extract amplitudes and phases
amplitudes=cell2mat(waveform(1))';
phases=cell2mat(waveform(2))';
initial_npoints=length(amplitudes);
initial_integral=simps(amplitudes);

% Do not resample by default
if nargin==1, npoints=initial_npoints; end

% Decide how to proceed
if npoints < initial_npoints
    
    % Decide the number of points to be truncated
    points_to_truncate=initial_npoints-npoints;
    if round(points_to_truncate/2)~=points_to_truncate/2
        error('read_wave: old and new point count must both be either odd or even.');
    end
    
    % Perform harmonic resampling
    complex_amplitudes=amplitudes.*cos(pi*phases/180)+1i*amplitudes.*sin(pi*phases/180);
    fourier_transform=fftshift(fft(complex_amplitudes));
    fourier_transform(1:points_to_truncate/2)=[]; fourier_transform(end-points_to_truncate/2+1:end)=[];
    complex_amplitudes=ifft(fftshift(fourier_transform));
    amplitudes=abs(complex_amplitudes);
    phases=90-(180/pi)*atan2(real(complex_amplitudes),imag(complex_amplitudes));
    integral=simps(amplitudes);
    
    % Print summary and diagnostics
    disp(['read_wave: waveform downsampled from ' num2str(initial_npoints) ' to ' num2str(npoints) ' points.']);
    disp(['read_wave: RF amplitude integral (Simpson quadrature) ' num2str(integral/initial_integral) ' of the original.' ]);
    
    % Bomb out if the integral changes by more than 5%
    if abs(integral/initial_integral-1)>0.05
        error('read_wave: integral check failure - the chosen downsampling level corrupts the waveform.');
    end
    
elseif npoints > initial_npoints
    
    % Could implement the upsampling here at some point
    error('read_wave: waveform upsampling not implemented.');
    
else
    
    % No resampling necessary.
    integral=initial_integral;
    
end

end

% I condemn Christianity... It is, to me, the greatest of all imaginable
% corruptions. [...] To breed in humans a self-contradiction, an art of
% self-pollution, a will to lie at any price, an aversion and contempt for
% all good and honest instincts! [...] The beyond as the will to deny all
% reality; the cross as the distinguishing mark of the most subterranean
% conspiracy ever heard of - against health, beauty, well-being,
% intellect... against life itself.
%
% I call Christianity the one great curse, the one great intrinsic
% depravity, the one great instinct of revenge, for which no means are
% venomous enough, or secret, subterranean and small enough - I call it 
% the one immortal blemish upon the human race...
%
% Friedrich Nietzsche

