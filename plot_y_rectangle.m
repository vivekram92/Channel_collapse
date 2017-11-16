n = 10;
p_hat = linspace(0,0.1,n);
w = 0.2;
h = 0.1;
v = 0.49;

for i = 1:n; 
%   y1(i) = triangle(p_hat(i)); 
   y_rectangle(i) = rectangle(p_hat(i),w/2,h/2,v); 
   i
%   y_rectangle(i) = rectangle(p_hat(i),0.1,0.05,0.5); 
%    Gauge = rectangle(p,w,h,v)/p; 
end

%y3_approx = (4/3)*eps*w/h;

A0 = w*h;

p = p_hat/(1-v^2);

A = w*h*(1 - p).*(1+v*(1+v)*p) - (4-(4*v*v))*sqrt(2)*p*w^2/3;

figure(1)
hold on
plot(p_hat,1-y_rectangle,'ro',p_hat,A/A0,'k--')
xlabel('Pressure, p')
