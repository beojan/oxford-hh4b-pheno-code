function y = CrystalBall2(x,xdata)

% x = vector of parameters to fit (alpha, n, xbar, sigma, Scale, n2, alpha2)
% n2/alpha2 may be other way round

%A = ( x(2) / abs(x(1)))^(x(2))*exp(-(abs(x(1))^2)/2);
%B = x(2)/abs(x(1)) - abs(x(1));
%N = 1/(sigma*(C + D));
%C = (x(2)/abs(x(1)))*(1/(x(2) - 1))*exp(-(abs(x(1))^2)/2);
%D = sqrt(pi/2)*(1 + erf(abs(x(1))/sqrt(2)));

y = x(5).*(heaviside((x(3)-x(1).*x(4)) - xdata).*(( x(2) ./ abs(x(1))).^(x(2)).*exp(-(abs(x(1)).^2)/2) .* (x(2)./abs(x(1)) - abs(x(1)) - (xdata - x(3))./x(4)).^-x(2)) + (heaviside(xdata - (x(3)-x(1).*x(4))).*heaviside((x(3)+x(7).*x(4)) - xdata).*(exp(-(xdata-x(3)).^2/(2.*x(4).^2)))) + (heaviside(xdata - (x(3)+x(7).*x(4))).*(( x(6) ./ abs(x(7))).^(x(6)).*exp(-(abs(x(7)).^2)/2) .* (x(6)./abs(x(7)) - abs(x(7)) - (x(3) - xdata)./x(4)).^-x(6))));
          
