clear;clf;

rng(111);

n=5;%n个随机插值节点

ui=pi*rand(n,1);
vi=2*pi*rand(n,1);
xi=sin(ui).*cos(vi);
yi=sin(ui).*sin(vi);
zi=cos(ui);
P=[xi';yi';zi'];

%绘制单位球面

t=linspace(0,pi,50);
p=linspace(0,2*pi,50);
[theta,phi]=meshgrid(t,p);
x=sin(theta).*sin(phi);
y=sin(theta).*cos(phi);
z=cos(theta);
re=[0 0 1];
colormap(re)
surf(x,y,z);
shading interp
axis equal;
daspect([1,1,1])
view(3);
axis tight
camlight
grid on
lighting gouraud
alpha(0.2);

%计算v，alpha

hold on
plot3(xi,yi,zi,'ro','LineWidth',2);
v=zeros(3,n-2);
for ii=2:n-1
    v_p=Log_sphere(P(:,ii),P(:,ii+1));
    v_p=v_p/norm(v_p);
    v_m=-Log_sphere(P(:,ii),P(:,ii-1));
    v_m=v_m/norm(v_m);
    v_temp=v_p+v_m;
    v(:,ii-1)=v_temp/norm(v_temp);
end

A=zeros(n-2,n-2);
c=zeros(n-2,1);
A(1,1:2)=[12*v(:,1)'*v(:,1) 3*v(:,2)'*v(:,1)];
c(1)=3*Log_sphere(P(:,2),P(:,3))'*v(:,1)-2*Log_sphere(P(:,2),P(:,1))'*v(:,1);
for ii=2:n-3
    A(ii,ii-1:ii+1)=[v(:,ii-1)'*v(:,ii) 4*v(:,ii)'*v(:,ii) v(:,ii+1)'*v(:,ii)];
    c(ii)=(Log_sphere(P(:,ii+1),P(:,ii+2))-Log_sphere(P(:,ii+1),P(:,ii)))'*v(:,ii);
end
A(n-2,n-3:n-2)=[3*v(:,n-3)'*v(:,n-2) 12*v(:,n-2)'*v(:,n-2)];
c(n-2)=3*Log_sphere(P(:,n-1),P(:,n-2))'*v(:,n-2)-2*Log_sphere(P(:,n-1),P(:,n))'*v(:,n-2);
al=A\c;


%计算控制点

B_p=zeros(3,n-2);B_m=B_p;
for ii=1:n-2
    if ii==1
        B_p(:,ii)=Exp_sphere(P(:,ii+1),al(ii)*v(:,ii));
        B_m(:,ii)=Exp_sphere(P(:,ii+1),-1.5*al(ii)*v(:,ii));
    else
        if ii==n-2
            B_p(:,ii)=Exp_sphere(P(:,ii+1),1.5*al(ii)*v(:,ii));
            B_m(:,ii)=Exp_sphere(P(:,ii+1),-al(ii)*v(:,ii));
        else
            B_p(:,ii)=Exp_sphere(P(:,ii+1),al(ii)*v(:,ii));
            B_m(:,ii)=Exp_sphere(P(:,ii+1),-al(ii)*v(:,ii));
        end
    end
end
B=[B_p,B_m];
plot3(B(1,:),B(2,:),B(3,:),'go','linewidth',2)


% 绘制曲线

t=linspace(0,1,25);
c=b2_sphere(P(:,1),B_m(:,1),P(:,2),t);
plot3(c(1,:),c(2,:),c(3,:),'r','LineWidth',2);
for ii=2:n-2
    c=b3_sphere(P(:,ii),B_p(:,ii-1),B_m(:,ii),P(:,ii+1),t);
    plot3(c(1,:),c(2,:),c(3,:),'r','LineWidth',2);
end
c=b2_sphere(P(:,n-1),B_p(:,n-2),P(:,n),t);
plot3(c(1,:),c(2,:),c(3,:),'r','LineWidth',2);

hold off
legend('单位球面','插值节点','控制点','贝塞尔曲线')



%单位球上的指数映射与对数映射

function p=Exp_sphere(p0,v)
s=norm(v);
p=p0*cos(s)+v/s*sin(s);
end

function v=Log_sphere(p0,p)
if norm(p0-p)<1e-10
    v=[0;0;0];
else
    theta=acos(p0'*p);
    v=-theta*cos(theta)*p0+theta*p;
    v=v/norm(v)*theta;
end
end

%单位球上的二次贝塞尔曲线

function c=b2_sphere(p0,p1,p2,t)
%p0,p1,p2为控制点(列向量)，t属于[0,1]
n=length(t);
c=zeros(3,n);
for ii=1:n
    theta01=acos(p0'*p1);
    c01=sin((1-t(ii))*theta01)/sin(theta01)*p0+sin(t(ii)*theta01)/sin(theta01)*p1;
    theta12=acos(p1'*p2);
    c12=sin((1-t(ii))*theta12)/sin(theta12)*p1+sin(t(ii)*theta12)/sin(theta12)*p2;
    theta02=acos(c01'*c12);
    c(:,ii)=sin((1-t(ii))*theta02)/sin(theta02)*c01+sin(t(ii)*theta02)/sin(theta02)*c12;
end
end

%单位球上的三次贝塞尔曲线

function c=b3_sphere(p0,p1,p2,p3,t)
%p0,p1,p2,p3为控制点(列向量)，t属于[0,1]
n=length(t);
c=zeros(3,n);
for ii=1:n
    theta01=acos(p0'*p1);
    c01=sin((1-t(ii))*theta01)/sin(theta01)*p0+sin(t(ii)*theta01)/sin(theta01)*p1;
    theta12=acos(p1'*p2);
    c12=sin((1-t(ii))*theta12)/sin(theta12)*p1+sin(t(ii)*theta12)/sin(theta12)*p2;
    theta23=acos(p2'*p3);
    c23=sin((1-t(ii))*theta23)/sin(theta23)*p2+sin(t(ii)*theta23)/sin(theta23)*p3;
    theta02=acos(c01'*c12);
    c02=sin((1-t(ii))*theta02)/sin(theta02)*c01+sin(t(ii)*theta02)/sin(theta02)*c12;
    theta13=acos(c12'*c23);
    c13=sin((1-t(ii))*theta13)/sin(theta13)*c12+sin(t(ii)*theta13)/sin(theta13)*c23;
    theta03=acos(c02'*c13);
    c(:,ii)=sin((1-t(ii))*theta03)/sin(theta03)*c02+sin(t(ii)*theta03)/sin(theta03)*c13;
end
end
