function [a,sigma] = burg(x, p)
    N = length(x);
    x =x(:);
         
    a = zeros(p+1,1);
    a(1) = 1;
    sigma = x'*x./N;
    f = x(2: end);
    b = x(1: end-1);
    for i = 1: p 
        k = -(2. * b' * f ) / (f' * f + b'*b);
        ftmp = f(2: end) + k .*b(2: end);
        b = b(1: end-1) + k' .* f(1:end-1);
        f = ftmp;	
        a(2: i+1) = a(2: i+1) + k .* conj(a(i: -1: 1));	
        sigma = (1 - k*k') * sigma;
    end 
   
end
