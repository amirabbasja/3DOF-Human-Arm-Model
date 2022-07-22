clc;
clear;

syms th1(t) th2(t) xd(t) T1(t) T2(t)
m1=1; m2 =1;md =1;g =1;c1 =1;c2 =1;cd =1;kd =1;l1 =1;l2 =1;l3 =1;lc1=1; lc2 =1;lc3=1; ld = 1;
%T = .5*m1*(lc1*dth1)^2+.5*m2*((l1*dth1)^2+(lc1*dth2)^2+2*l1*lc2*dth1*dth2*cos(th1-th2))+.5*md*(dl*dth1+dxd)^2;
%V = m1*g*lc1*sin(th1)+m2*g*(l1*sin(th1)+lc2*sin(th2))+md*(dl*sin(th1)+xd*cos(th1))+.5*kd*xd^2;

eq1 = diff(th1, 2)*(m1*lc1^2+m2*l1^2+md*ld^2) + m2*diff(th2,2)*l1*lc2*cos(th2-th1) -l1*lc2*sin(th2-th1)*diff(th2,1)^2 ...
    + md*ld*diff(xd,2)+m1*g*lc1*cos(th1) + m2*g*l1*cos(th1) + md*(ld*cos(th1)-xd*sin(th1)) == sin(t) - c1*diff(th1,1);

eq2 = m2*lc2^2*diff(th2,2) + m2*l1*lc2*diff(th1,2)*cos(th2-th1)+m2*l1*lc2*diff(th1,1)^2*sin(th2-th1) ...
    +m2*g*lc2*cos(th2) == cos(t) - c2*diff(th2,1);

eq3 = md*(diff(xd,2)+ld*diff(th1,2)) + md*cos(th1) + kd*xd == 0 - cd*diff(xd,1);

[V,S] = odeToVectorField(eq1, eq2, eq3);
M = matlabFunction(V,'vars', {'t','Y'});

  xd_0 = 0;
 Dxd_0 = 0;
 th1_0 = 0;
Dth1_0 = 0;
 th2_0 = 0;
Dth2_0 = 0;

sol = ode45(M,[0 40],[xd_0 Dxd_0 th1_0 Dth1_0 th2_0 Dth2_0]);

x = linspace(0,40,100);
y = deval(sol,x,1);
plot(x,y)

