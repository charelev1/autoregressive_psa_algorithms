function [a, sigma] = burg_fast_fft(x, p)
    N = length(x);
    x = x(:);
   	a = zeros(p+1, 1); % prediction error coeffs
    c = zeros(p+1, 1); % auto correlation coeffs
    
    % FFT-based autocorrelation calculation
    xf = [x(:); zeros(2^nextpow2(N) - N ,1)]; % zero pad and FFT
    fx1 = fft(xf); 
    fxsqr = fx1.*conj(fx1);
    c = ifft(fxsqr);
    c = c(1: p+1);

    %for j = 1: p+1 
    %    c(j) = x(1: N-j+1)' *  x(j : N);
    %end
    
    a(1) = 1;
    g = [2 * c(1) - abs(x(1)) ^ 2 - abs(x(N)) ^ 2; 2 * c(2)];
    r = 2 * c(2);
    sigma = c(1)./N;
    
    %% intermediate steps
    for i = 1: p-1
        k = - (a(1: i)' * g(end: -1: 2)) / (a(1: i)' * g(1:end -1));
        a(2: i + 1) = a(2: i + 1) + k .* conj(a(i: -1: 1));	
        r = [2 * c(i+2); r - x(1:i) * x(i+1) - x(N: -1: N-i+1) * x(N-i)];
        DR  = - x(i+1: -1: 1) * x(i+1: -1: 1)' - x(N-i: N) * x(N-i: N)'; 
        g = [g + k * g(end: -1: 1) + DR * a(1: i+1); r' * a(1: i+1)];
		sigma = (1 - k * k') * sigma;
    end 
    
    %% final step
    k = - (a(1: p)' * g(end: -1: 2)) / (a(1: p)' * g(1: end -1));
    a(2: p+1) = a(2: p+1) + k .* conj(a(p: -1: 1));	
	sigma = (1 - k * k') * sigma;
end
