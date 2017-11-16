clear
clc

p = 0.01;

w = 0.1;
h = [0.01:0.01:0.09 0.0999];
%h = [0.01:0.01:0.1];
v = 0.49;

eps = p/(1-v^2);
n = length(h);
for i = 1:n; 
    Gauge(i) = rectangle(p,w,h(i),v)/eps;
    i
end

aspect = linspace(0.06,1.1,1000);

%Gauge_approx = (4/3)./aspect + (1 - (1-eps)*(1+v*eps))/eps;
Gauge_approx = 1 - v - v^2 + (4 - (4*v*v))*sqrt(2)./(3*aspect);

figure(2);  hold on
plot(h/w,Gauge,'rs',aspect,Gauge_approx,'k--')
% figure(2);  hold on
% plot(aspect,Gauge_approx,'k--')