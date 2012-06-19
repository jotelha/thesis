% Iterative thresholding baseline level estimator.
%
% ilya.kuprov@oerc.ox.ac.uk

function baseline_level=baseline(spectrum,threshold,n_iter)

spectrum_real=real(spectrum(:));
spectrum_imag=imag(spectrum(:));

baseline_real=true(size(spectrum_real));
baseline_imag=true(size(spectrum_imag));

for n=1:n_iter
    mean_real=mean(spectrum_real(baseline_real));
    mean_imag=mean(spectrum_imag(baseline_imag));
    stdev_real=std(spectrum_real(baseline_real));
    stdev_imag=std(spectrum_imag(baseline_imag));
    baseline_real=(abs(spectrum_real-mean_real)<threshold*stdev_real);
    baseline_imag=(abs(spectrum_imag-mean_imag)<threshold*stdev_imag);
end

baseline_level=mean(spectrum_real(baseline_real))+1i*mean(spectrum_imag(baseline_imag));

end

% ...and it seemed so simple, so unastonishing, that the thing she felt was
% like a blessing pronounced upon the universe by means of three words: but
% of course.
%
% Ayn Rand, "Atlas Shrugged"

