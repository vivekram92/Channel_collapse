function y = rectangle(P,W,H,nu);
%%% function p0 = rectangle(0.07,0.1,0.5,0.4);

model = createpde(2);

L1 = 10;
L2 = 1;
w = W;
h = H;

v = nu;

p0 = P;

R1 = [3,4,-L1/2,L1/2,L1/2,-L1/2,-L2/2,-L2/2,L2/2,L2/2]';
R2 = [3,4,-w/2,w/2,w/2,-w/2,-h/2,-h/2,h/2,h/2]';
R2 = [R2;zeros(length(R1)-length(R2),1)];

gm = [R1,R2];
ns = (char('R1','R2'))';
sf = 'R1-R2';
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
% c = char('1','0','0',chi,'0','0',chi_nu,'0','0','0',chi_nu,'0',chi,'0','0','1');
c = char('1','0','0',chi,'0','0',chi_nu,'0','0','0',chi_nu,'0',chi,'0','0','1');
% c = char('1','0','0','0','0','0.5','0.5','0','0','0.5','0.5','0','0','0','0','1');

a = char('0','0');
f = char('0','0');

Q = zeros(2,2);
G1 = [0;-p0];
G2 = [0;p0];
G3 = [0;0];
applyBoundaryCondition(model,'Edge',[2],'q',Q,'g',G1,'Vectorized','on');
applyBoundaryCondition(model,'Edge',[7],'q',Q,'g',G2,'Vectorized','on');
applyBoundaryCondition(model,'Edge',[1,3,4,5,6,8],'q',Q,'g',G3,'Vectorized','on');
%applyBoundaryCondition(model,'Edge',[6],'u',[0,0],'Vectorized','on');

% mesh = generateMesh(model);
% u = pdenonlin(model,c,a,f);
mesh = generateMesh(model,'Hmax',0.7e-2);
u = pdenonlin(model,c,a,f,'Tol',1e-8);


N = length(u)/2;

n = 100;
% X = linspace(-L1/2,L1/2,n);
% Y = linspace(-L2/2,L2/2,n);

[p,e,t] = meshToPet(mesh);
% u1 = tri2grid(p,t,u(1:N),X,Y);
% u2 = tri2grid(p,t,u((N+1):2*N),X,Y);
% 
% XL = -L1/2 + u1(1,n/2);
% YL = u2(1,n/2);
% XR = L1/2 + u1(n,n/2);
% YR = u2(n,n/2);
% 
% u0_X = (XL + XR)/2;
% u0_Y = (YL + YR)/2;
% 
% figure(4); hold on
% plot(X+u1(n,:)-u0_X*ones(1,n),(L2/2)*ones(1,n)+u2(n,:)-u0_Y*ones(1,n),'k-')
% plot(X+u1(1,:)-u0_X*ones(1,n),-(L2/2)*ones(1,n)+u2(1,:)-u0_Y*ones(1,n),'k-')
% plot((L1/2)*ones(n,1)+u1(:,n)-u0_X*ones(n,1),Y'+u2(:,n)-u0_Y*ones(n,1),'k-')
% plot(-(L1/2)*ones(n,1)+u1(:,1)-u0_X*ones(n,1),Y'+u2(:,1)-u0_Y*ones(n,1),'k-')
% 
% plot(X,(L2/2)*ones(1,n),'k--')
% plot(X,-(L2/2)*ones(1,n),'k--')
% plot((L1/2)*ones(n,1),Y','k--')
% plot(-(L1/2)*ones(n,1),Y','k--')


% x = linspace(-w/2,w/2,n);
% y = linspace(-h/2,h/2,n);
% [p,e,t] = meshToPet(mesh);
% u1 = tri2grid(p,t,u(1:N),x,y);
% u2 = tri2grid(p,t,u((N+1):2*N),x,y);
% 
% figure(4); hold on
% plot(x+u1(n,:)-u0_X*ones(1,n),(h/2)*ones(1,n)+u2(n,:)-u0_Y*ones(1,n),'k-')
% plot(x+u1(1,:)-u0_X*ones(1,n),-(h/2)*ones(1,n)+u2(1,:)-u0_Y*ones(1,n),'k-')
% plot((w/2)*ones(n,1)+u1(:,n)-u0_X*ones(n,1),y'+u2(:,n)-u0_Y*ones(n,1),'k-')
% plot(-(w/2)*ones(n,1)+u1(:,1)-u0_X*ones(n,1),y'+u2(:,1)-u0_Y*ones(n,1),'k-')
% 
% plot(x,(h/2)*ones(1,n),'k--')
% plot(x,-(h/2)*ones(1,n),'k--')
% plot((w/2)*ones(n,1),y','k--')
% plot(-(w/2)*ones(n,1),y','k--')
% %axis equal
% axis([-0.6*L2 0.6*L2 -0.6*L2 0.6*L2])
% axis square

% figure(2)
% surf(X,Y,u1,'FaceAlpha',0.5,'EdgeColor','none')
% xlabel('X_1 (mm)')
% ylabel('X_2 (mm)')
% zlabel('u_1 (mm)')
% colormap bone
% axis equal
% 
% figure(3)
% surf(X,Y,u2,'FaceAlpha',0.5,'EdgeColor','none')
% xlabel('X_1 (mm)')
% ylabel('X_2 (mm)')
% zlabel('u_2 (mm)')
% colormap bone
% axis equal

F1 = pdeInterpolant(p,t,u(1:N));
F2 = pdeInterpolant(p,t,u((N+1):2*N));

x = [linspace(-w/2,-w/2,n) linspace(-w/2,w/2,n) linspace(w/2,w/2,n) linspace(w/2,-w/2,n)];
dx = diff(x); dx = [dx(1) dx];
y = [linspace(-h/2,h/2,n) linspace(h/2,h/2,n) linspace(h/2,-h/2,n) linspace(-h/2,-h/2,n)];
dy = diff(y); dy = [dy(1) dy];

for i = 1:4*n
    uX = evaluate(F1,x(i),y(i));
    uY = evaluate(F2,x(i),y(i));    
    X_t(i) = x(i) + uX;
    Y_t(i) = y(i) + uY;
end

A0 = -0.5*sum(x.*dy - y.*dx);

dX_t = diff(X_t); dX_t = [dX_t(1) dX_t];
dY_t = diff(Y_t); dY_t = [dY_t(1) dY_t];
Area = -0.5*sum(X_t.*dY_t - Y_t.*dX_t);

%eps = p0/(1-v^2);
%Ap = A0*(1 - eps)*(1+v*eps) - (4/3)*eps*w^2;
%Ap = A0 - (4/3)*eps*w^2;
%[A/A0 Ap/A0]
y = 1 - Area/A0;
