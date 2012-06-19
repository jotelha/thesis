% Contour plotting utility with non-linear adaptive contour spacing. The
% function is useful for NMR data where small cross-peaks must be
% adequately contoured next to large diagonal peaks. The following function
% is used to define contour level positions:
%
% positive_contours=xmax*(1-2*delta)*((0:(ncont-1))/(ncont-1)).^k+xmax*delta
%
% where:
%
%       * xmax and xmin are calculated from the spectrum;
%
%       * delta is the minimum and maximum elevation (as a fraction of total
%         intensity) of the contours above the baseline. A reasonable value for 
%         most 2D spectra is [0.02 0.2];
%       
%       * ncont is the number of contours;
%
%       * k controls the curvature of the contour spacing function: k=1 
%         corresponds to linear spacing and k>1 bends the spacing curve to
%         increase the sampling density near the baseline. A reasonable
%         value is 2;
%
%       * ncol is a number of colors in the colormap (around 256 is fine);
%
%       * m is the curvature of the colormap: m=1 corresponds to a linear
%         color ramp into the red for positive contours and into the blue
%         for negative contours. A reasonable value for high-contrast
%         plotting is 6.
%
%       * signs can be set to 'positive', 'neagtive' or 'both' - this will
%         cause the corresponding contours to be plotted.
%
% ilya.kuprov@oerc.ox.ac.uk
% matthew.krzystyniak@oerc.ox.ac.uk

function contour_plot(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)

% Set the defaults
if nargin==3
    ncont=20; delta=[0.02 1.0]; k=2; ncol=256; m=6; signs='both';
elseif nargin==8
    signs='both';
end

% Inspect the spectrum
if nnz(imag(spectrum))>0
    spectrum=real(spectrum);
    report(spin_system,'contour_plot: complex spectrum received, plotting the real part...');
end

% Determine the data extents
xmax=max(spectrum(:)); xmin=min(spectrum(:)); 

% Compute the contour levels
if (xmax>0)&&(strcmp(signs,'positive')||strcmp(signs,'both'))
    positive_contours=delta(2)*xmax*linspace(0,1,ncont).^k+xmax*delta(1);
else
    positive_contours=[];
end
if (xmin<0)&&(strcmp(signs,'negative')||strcmp(signs,'both'))
    negative_contours=delta(2)*xmin*linspace(0,1,ncont).^k+xmin*delta(1);
else
    negative_contours=[];
end
contours=[negative_contours(end:-1:1) positive_contours];

% Get the axes
[a_f1,a_f2]=axis_2d(spin_system,parameters);

% Plot the spectrum
contour(a_f1,a_f2,spectrum,contours);

% Invert the axes
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');

% Label the axes
xlabel(parameters.axis_units);
ylabel(parameters.axis_units);

% Colour the contours
if any(positive_contours)&&any(negative_contours)
    plot_range=max(abs(positive_contours))+max(abs(negative_contours));
    nred=ceil(ncol*max(abs(positive_contours))/plot_range);
    nblue=ceil(ncol*max(abs(negative_contours))/plot_range);
elseif any(positive_contours)
    plot_range=max(abs(positive_contours));
    nred=ceil(ncol*max(abs(positive_contours))/plot_range);
    nblue=0;
elseif any(negative_contours)
    plot_range=max(abs(negative_contours));
    nblue=ceil(ncol*max(abs(negative_contours))/plot_range);
    nred=0;
else
    error('contour_plot: spectrum contouring produced no contours.');
end

colors=0.9*(1-[linspace(1,0,nblue)' linspace(1,0,nblue)' linspace(0,0,nblue)';
               linspace(0,0,nred)'  linspace(0,1,nred)'  linspace(0,1,nred)']).^m;
colormap(colors); colorbar; drawnow;

end

% After all, every murderer when he kills runs the risk of the most
% dreadful of deaths, whereas those who kill him risk nothing except
% promotion.
%
% Albert Camus 

