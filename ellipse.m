function y = ellipse(p,A,B,v);
%%% function p0 = ellipse(0.07);

model = createpde(2);

L1 = 10;
L2 = 1;

p0 = p;

R1 = [3,4,-L1/2,L1/2,L1/2,-L1/2,-L2/2,-L2/2,L2/2,L2/2]';
C1 = [4,0,0,A,B,0]';
C1 = [C1;zeros(length(R1)-length(C1),1)];

gm = [R1,C1];
ns = (char('R1','C1'))';
sf = 'R1-C1';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);

% figure(1)
% pdegplot(g,'EdgeLabels','on')
% axis equal

% nu = sprintf('%d',v);
% chi = sprintf('(1-%d)/4',v);
nu = sprintf('(%d)/(1-%d)',v,v);
chi = sprintf('(1-(2*(%d)))/(2*(1-(%d)))',v,v);
chi_nu = sprintf('((1-(2*(%d)))/(2*(1-(%d)))) + ((%d)/(1-%d))',v,v,v,v);

% c = char('1','0','0',chi,'0',nu,chi,'0','0',chi,nu,'0',chi,'0','0','1');
% c = char('1','0','0',chi,'0',chi,nu,'0','0',nu,chi,'0',chi,'0','0','1');
c = char('1','0','0',chi,'0','0',chi_nu,'0','0','0',chi_nu,'0',chi,'0','0','1');
% c = char('1','0','0','0','0','0.5','0.5','0','0','0.5','0.5','0','0','0','0','1');

a = char('0','0');
f = char('0','0');

Q = zeros(2,2);
G1 = [0;-p0];
G2 = [0;p0];
G3 = [0;0];
applyBoundaryCondition(model,'Edge',[2],'q',Q,'g',G1,'Vectorized','on');
applyBoundaryCondition(model,'Edge',[4],'q',Q,'g',G2,'Vectorized','on');
applyBoundaryCondition(model,'Edge',[1,5,6,7,8],'q',Q,'g',G3,'Vectorized','on');
applyBoundaryCondition(model,'Edge',[3],'u',[0,0],'Vectorized','on');

% mesh = generateMesh(model);
% u = pdenonlin(model,c,a,f);
mesh = generateMesh(model,'Hmax',0.007);
u = pdenonlin(model,c,a,f,'Tol',1e-8);
[p,e,t] = meshToPet(mesh);

N = length(u)/2;
n = 500;

F1 = pdeInterpolant(p,t,u(1:N));
F2 = pdeInterpolant(p,t,u((N+1):2*N));
theta = linspace(0,2*pi,n);
for i = 1:n
    x(i) = A*sin(theta(i));
    y(i) = B*cos(theta(i));
    uX = evaluate(F1,x(i),y(i));
    uY = evaluate(F2,x(i),y(i));

    X_R(i) = x(i) + uX;
    Y_R(i) = y(i) + uY;

end

dx = diff(x); dx = [dx(1) dx];
dy = diff(y); dy = [dy(1) dy];
A0 = -0.5*sum(x.*dy - y.*dx);

dX_R = diff(X_R); dX_R = [dX_R(1) dX_R];
dY_R = diff(Y_R); dY_R = [dY_R(1) dY_R];
Area = -0.5*sum(X_R.*dY_R - Y_R.*dX_R);
[pi*A*B/A0]
y = 1 - Area/A0;


