clear
clc

p = 0.01;

A = 0.1;
B = [0.01:0.01:0.1];
v = 0.49;

eps = p/(1-v^2);
n = length(B);
for i = 1:n; 
    Gauge(i) = ellipse(p,A,B(i),v)/eps;
    i
end

figure(1); hold on
plot(B/A,Gauge,'ro')


aspect = linspace(0.08,1.1,1000);
B = A*aspect;
A0 = pi*A*B;

% k = (3-v)/(1+v); 
k = (3-(4*v)); 
chi = (k+1)/8;
S0 = (1 - B/A)./(chi*(3 + B/A));
S = 2*(1+v)*p/(1-v^2) + S0;
R0 = A./(1 + S0*chi);

Ap = pi*(R0.^2).*(1 + chi*S).*(1-3*chi*S);

%Gauge_approx = (1 - Ap./A0)/eps;
Gauge_approx = (9 + aspect.^2)./(4*aspect);

figure(1); hold on
plot(aspect,Gauge_approx,'k--')
xlabel('Aspect Ratio, \alpha')
ylabel('Gauge Factor, G_e')


