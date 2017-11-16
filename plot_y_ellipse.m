n = 10;
p_hat = linspace(0,0.1,n);
w = 0.2;
h = 0.1;
v = 0.49;

for i = 1:n; 
%   y1(i) = triangle(p_hat(i)); 
   y_ellipse(i) = ellipse(p_hat(i),w/2,h/2,v); 
   i
%   y_rectangle(i) = rectangle(p_hat(i),0.1,0.05,0.5); 
%    Gauge = rectangle(p,w,h,v)/p; 
end

%y3_approx = (4/3)*eps*w/h;

A = w/2;
B = h/2;
A0 = pi*A*B;

% k = (3-v)/(1+v); 
k = (3-(4*v)); 
chi = (k+1)/8;
p0 = (1 - B/A)./(2*chi*(1+v)*(3 + B/A));
R0 = A/(1 + 2*(1+v)*chi*p0);

p = p_hat/(1-v^2);

q = chi*(1+v)*(p + p0);
Ap = pi*(R0.^2).*(1 + 2*q).*(1-6*q);

%S0 = (1 - B/A)./(chi*(3 + B/A));
% S = 2*(1+v)*p_hat/(1-v^2) + S0;
% R0 = A./(1 + S0*chi);
%Ap = pi*(R0.^2).*(1 + chi*S).*(1-3*chi*S);

figure(2)
plot(p_hat,1-y_ellipse,'ro',p_hat,Ap/A0,'k--')
