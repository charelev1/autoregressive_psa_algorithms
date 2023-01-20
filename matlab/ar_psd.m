function PSD = ar_psd(sigma, a, p, w)
% Estimate the PSD of a p-order AR model at the w normalized frequencies
    a = a(:);
	PSD = zeros(length(w),1); % output PSD

    %for each freq point compute the PSD
    for i = 1: length(w)
        ee = exp(-1i*(0:p)*w(i)); %exponent part (complex unity phazor)
        dem = dot(ee, a); % dot product a summation
        PSD(i) = 1/( dem*dem'); % squared denumenator
    end

    PSD = PSD*sigma;
	PSD = PSD(1: end/2 + 1);
end
