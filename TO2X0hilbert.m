function [y] = TO2X0(v, x, N1, N2, G)
gma_ho = v(1);
Gma_in = G;%v(9);
A1 = v(2);
x01 = v(3); 
y01 = v(4);
x1 = x(N1+1:N1+N2);

A2 = v(5);
x02 = v(6);
x2 = x(1:N1);

phi = v(7);
om = v(8);
phi2 = v(9);

%N = length(x1);

% Definition of Voigt from Wikipedia: http://en.wikipedia.org/wiki/Voigt_profile
y(N1+1:N1+N2) = real((y01 + A1*real(cerf(((x1-x01)+1i*gma_ho)./(Gma_in*sqrt(2))))/(Gma_in*sqrt(2*pi)))*exp(-1i*phi));
y(2*N1+N2+1:2*N1+2*N2) = imag((y01 + A1*real(cerf(((x1-x01)+1i*gma_ho)./(Gma_in*sqrt(2))))/(Gma_in*sqrt(2*pi)))*exp(-1i*phi)); %fits imaginary part

y(1:N1) = real(hilbert(A2*abs(exp((gma_ho - 1i*(x2-x02)).^2./(2*Gma_in^2)).*(1-erfz((gma_ho - 1i*(x2-x02))./(sqrt(2)*Gma_in)))...
    ./(Gma_in.*(gma_ho - 1i*(x2-x02))))).*exp(-1i*(phi+phi2+om*(x2-x02))));
y(N1+N2+1:2*N1+N2) = imag(hilbert(A2*abs(exp((gma_ho - 1i*(x2-x02)).^2./(2*Gma_in^2)).*(1-erfz((gma_ho - 1i*(x2-x02))./(sqrt(2)*Gma_in)))...
    ./(Gma_in.*(gma_ho - 1i*(x2-x02))))).*exp(-1i*(phi+2.5-0.58*(x2-x02)))); %hilbert transform breaks it into real and imaginary parts. 
