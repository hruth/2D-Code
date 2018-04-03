function [y] = TO2X0(v, x, N1, N2, G)
gma_ho = v(1); % guess for homo linewidth
Gma_in = v(7);%v(7); %guess for inhomogenous linewidth
A1 = v(2); % an amplitude guess (currently pkRatio2 should be pkRatio1)
x01 = v(3); % guess for the center of the inhomo .. should change to linewidths.LorentzianCenter
y01 = v(4); % 0 .. an offset guess. Should we really be fitting this?
x1 = x(N1+1:N1+N2);

A2 = v(5); % an amplitude guess from lorentzian
x02 = v(6); %center of Lorentzian
x2 = x(1:N1);

%N = length(x1);

% Definition of Voigt from Wikipedia: http://en.wikipedia.org/wiki/Voigt_profile
y(N1+1:N1+N2) = y01 + A1*real(cerf(((x1-x01)+i*gma_ho)./(Gma_in*sqrt(2))))/(Gma_in*sqrt(2*pi));

y(1:N1) = A2*abs(exp((gma_ho - i*(x2-x02)).^2./(2*Gma_in^2)).*(1-erfz((gma_ho - i*(x2-x02))./(sqrt(2)*Gma_in)))...
    ./(Gma_in.*(gma_ho - i*(x2-x02))));
