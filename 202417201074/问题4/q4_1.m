p = 1.7; % 螺距
a = p / (2 * pi); % 螺线系数
d1 = 2.86; % 龙头两把手间距
d2 = 1.65; % 除龙头板外其他板两把手间距
v0 = 1; % 龙头运动速度
R=4.5; % 掉头区域半径

% 绘制中心对称螺线
theta=5*2*pi:-0.01:0*pi;
r=a*theta;
x=r.*cos(theta);
y=r.*sin(theta);
figure(1)
hold on
set(gcf,'Position',[200 200 600 600]);
title('盘入、盘出螺线及掉头圆弧')
plot(x,y,'-','Color',[0.5, 0.5, 0.5],'DisplayName','盘入螺线')
axis equal
grid on
xlabel('x轴/m')
ylabel('y轴/m')

 % 中心对称图形，方位角相差180°
theta=theta-pi;
ra=a*(theta+pi);
x2=ra.*cos(theta);
y2=ra.*sin(theta);
plot(x2,y2,'m','Color','black','DisplayName','盘出螺线')
% 绘制掉头区域
x_diao=R*cos(theta);
y_diao=R*sin(theta);
plot(x_diao,y_diao,'Color','red','LineStyle','--','LineWidth',2,'DisplayName','掉头区域')

% 求两段圆弧的位置、切入切出点
theta_ru=R/a;
theta_chu=R/a-pi; % 另一个螺线
slope=(a*sin(theta_ru)+R*cos(theta_ru))/(a*cos(theta_ru)-R*sin(theta_ru));
theta_max_a=atan(-1/slope)+pi;
theta_equal=atan(tan(theta_ru))+pi-theta_max_a;
rC1_C2=R/cos(theta_equal);
r2=rC1_C2/3;
r1=r2*2;
phi=2*theta_equal;
SC1=r1*(pi-phi); SC2=r2*(pi-phi);
theta_min_a=theta_max_a-SC1/r1;
theta_min_2=theta_min_a-pi;
theta_max_2=theta_min_2+SC2/r2;
% 两段弧圆心坐标
x_C1=R*cos(theta_ru)+r1*cos(theta_max_a-pi);
y_C1=R*sin(theta_ru)+r1*sin(theta_max_a-pi);
x_C2=R*cos(theta_chu)-r2*cos(theta_max_2);
y_C2=R*sin(theta_chu)-r2*sin(theta_max_2);

% 继续完善图，添加绘制两段弧
hold on
plot(x_C1+r1*cos(linspace(theta_min_a,theta_max_a,50)),y_C1+r1*sin(linspace(theta_min_a,theta_max_a,50)),'Color','magenta','LineWidth',2)
plot(x_C1,y_C1,'ko','MarkerSize',5,'MarkerFaceColor','magenta')
plot(x_C2+r2*cos(linspace(theta_min_2,theta_max_2,50)),y_C2+r2*sin(linspace(theta_min_2,theta_max_2,50)),'Color','magenta','LineWidth',2)
plot(x_C2,y_C2,'ko','MarkerSize',5,'MarkerFaceColor','magenta')
axis equal
set(gca, 'looseInset', [0 0 0 0]);

% 盘入曲线上头节点的位置求解
my_theta=@(t,theta)1./(a*sqrt(1+theta.^2));
theta0=theta_ru;
dt=0.1; % 时间步长
t_span=0:dt:100;
[ttl,theta]=ode45(my_theta,t_span,theta0); % 龙格库塔法求解
X1=a*theta.*cos(theta);
Y1=a*theta.*sin(theta);
tt_ru=ttl(end:-1:1);
X=zeros(224,200/dt+1);
Y=zeros(224,200/dt+1);
Ting=zeros(224,200/dt+1);
X(1,1:length(X1))=X1(end:-1:1);
Y(1,1:length(Y1))=Y1(end:-1:1);
Ting(1,1:length(theta))=theta(end:-1:1);
% 圆弧上切点的位置信息
td_c1=dt:dt:SC1;
theta_C1=-td_c1/r1+theta_max_a;
Ting(1,length(theta)+(1:length(td_c1)))=theta_C1;
X(1,length(X1)+(1:length(td_c1)))=r1*cos(theta_C1)+x_C1;
Y(1,length(Y1)+(1:length(td_c1)))=r1*sin(theta_C1)+y_C1;
td_c2=td_c1(end)+dt:dt:SC1+SC2;
theta_C2=(td_c2-SC1)/r2+theta_min_a-pi;
Ting(1,length(theta)+length(theta_C1)+(1:length(td_c2)))=theta_C2;
X(1,length(X1)+length(theta_C1)+(1:length(td_c2)))=r2*cos(theta_C2)+x_C2;
Y(1,length(Y1)+length(theta_C1)+(1:length(td_c2)))=r2*sin(theta_C2)+y_C2;
my_theta=@(t,theta)1./(a*sqrt(1+(theta+pi).^2));
theta0=theta_chu;
t_span=td_c2(end)+dt:dt:100;
[ttl,theta2]=ode45(my_theta,t_span,theta0); % 龙格库塔法求解
X2=a*(theta2+pi).*cos(theta2);
Y2=a*(theta2+pi).*sin(theta2);
Ting(1,length(theta)+length(theta_C1)+length(td_c2)+1:end)=theta2;
X(1,length(theta)+length(theta_C1)+length(td_c2)+1:end)=X2;
Y(1,length(theta)+length(theta_C1)+length(td_c2)+1:end)=Y2;
set(gcf,'Position',[300 300 600 600])
waitfor(gcf);

% 绘制第二张图
figure(2)
for i=1:3:length(Ting(1,:))
    title('头把手中心的轨迹(-100到100s)')
    xlabel('x轴/m')
    ylabel('y轴/m')
    plot(X(1,i),Y(1,i),'Marker','o','MarkerSize',2,'MarkerFaceColor','r','Color','cyan')
    hold on
    axis equal
    axis([-10 10 -10 10])
    grid on
    drawnow
end
set(gca, 'looseInset', [0 0 0 0]);
waitfor(gcf);

% 求各段龙身的位置
N=223;
t_total=-100:dt:100;
for j=1:length(t_total)
    if t_total(j)<=0
        for i=2:N+1
            d=d1*(i<=2)+d2*(i>2);
            theta_xy=solve_theta1(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d);
            Ting(i,j)=theta_xy;
            X(i,j)=a*theta_xy*cos(theta_xy);
            Y(i,j)=a*theta_xy*sin(theta_xy);
        end
    elseif t_total(j)>0 && t_total(j)<=SC1
        f=2;
        for i=2:N+1
            d=d1*(i<=2)+d2*(i>2);
            if f==2
                [xi,yi,arf,f]=solve_point_2_1(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d,r1,x_C1,y_C1,theta_max_a);
                Ting(i,j)=arf;
                X(i,j)=xi;
                Y(i,j)=yi;
            else
                theta_xy=solve_theta1(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d);
                Ting(i,j)=theta_xy;
                X(i,j)=a*theta_xy*cos(theta_xy);
                Y(i,j)=a*theta_xy*sin(theta_xy);
            end
        end
    elseif t_total(j)>SC1 && t_total(j)<=SC1+SC2
        f=3;
        for i=2:N+1
            d=d1*(i<=2)+d2*(i>2);
            if f==3
                [xi,yi,arf,f]=solve_point_3_2(X(i-1,j),Y(i-1,j),Ting(i-1,j),d,r1,x_C1,y_C1,r2,x_C2,y_C2,theta_min_2);
                Ting(i,j)=arf;
                X(i,j)=xi;
                Y(i,j)=yi;
            elseif f==2
                [xi,yi,arf,f]=solve_point_2_1(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d,r1,x_C1,y_C1,theta_max_a);
                Ting(i,j)=arf;
                X(i,j)=xi;
                Y(i,j)=yi;
            else
                theta_xy=solve_theta1(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d);
                Ting(i,j)=theta_xy;
                X(i,j)=a*theta_xy*cos(theta_xy);
                Y(i,j)=a*theta_xy*sin(theta_xy);
            end
        end
    else
        f=4;
        for i=2:N+1
            d=d1*(i<=2)+d2*(i>2);
            if f==4
                [xi,yi,arf,f]=solve_point_4_3(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d,r2,x_C2,y_C2,theta_max_2);
                Ting(i,j)=arf;
                X(i,j)=xi;
                Y(i,j)=yi;
            elseif f==3
                [xi,yi,arf,f]=solve_point_3_2(X(i-1,j),Y(i-1,j),Ting(i-1,j),d,r1,x_C1,y_C1,r2,x_C2,y_C2,theta_min_2);
                Ting(i,j)=arf;
                X(i,j)=xi;
                Y(i,j)=yi;
            elseif f==2
                [xi,yi,arf,f]=solve_point_2_1(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d,r1,x_C1,y_C1,theta_max_a);
                Ting(i,j)=arf;
                X(i,j)=xi;
                Y(i,j)=yi;
            else
                theta_xy=solve_theta1(p,X(i-1,j),Y(i-1,j),Ting(i-1,j),d);
                Ting(i,j)=theta_xy;
                X(i,j)=a*theta_xy*cos(theta_xy);
                Y(i,j)=a*theta_xy*sin(theta_xy);
            end
        end
    end
end

% 调头过程可视化
figure(3)
clf;
set(gcf,'Position',[200 200 600 600])
for j=1:size(X,2)
    plot(X(:,j),Y(:,j),'k-',Marker','*','MarkerSize',4,'MarkerFaceColor','cyan')
    title('100s时板凳龙的轨迹')
    axis equal
    grid on
    xlabel('x轴/m')
    ylabel('y轴/m')
    axis([-25 25 -25 25])
    drawnow
end
set(gca, 'looseInset', [0 0 0 0]);

% 保存至相应位置的excel文件
xlswrite('X.xlsx',X)
xlswrite('Y.xlsx',Y)

function theta=solve_theta1(p,x,y,theta,d)
a=p/2/pi;
fun=@(theta)(a*theta.*cos(theta)-x).^2+(a*theta.*sin(theta)-y).^2-d^2;
q=0.01;
options = optimoptions('fsolve','Display','off');
theta=fsolve(fun,theta+q,options);
while theta<=theta || abs(a*theta-a*theta)>a/2
    q=q+0.1;
    theta=fsolve(fun,theta+q,options);
end
end

function theta=solve_theta2(p,x,y,theta,d)
a=p/2/pi;
fun=@(theta)(a*(theta+pi).*cos(theta)-x).^2+(a*(theta+pi).*sin(theta)-y).^2-d^2;
q=-0.1;
options = optimoptions('fsolve','Display','off');
theta=fsolve(fun,theta+q,options);
while theta>=theta || abs(a*theta-a*theta)>p/2
    q=q-0.1;
    theta=fsolve(fun,theta+q,options);
end
end

function [x,y,theta,flag]=solve_point_2_1(p,x,y,theta,d,r1,xc,yc,theta_m)
a=p/2/pi;
delta_theta=2*asin(d/2/r1);
if delta_theta<=theta_m-theta
    flag=2;
    theta=theta+delta_theta;
    x=xc+r1*cos(theta);
    y=yc+r1*sin(theta);
else
    theta=solve_theta1(p,x,y,4.5/a,d);
    flag=1;
    x=a*theta*cos(theta);
    y=a*theta*sin(theta);
end
end

function [x,y,theta,flag]=solve_point_3_2(x,y,theta,d,r1,x1,y1,r2,x2,y2,theta_m)
delta_theta=2*asin(d/2/r2);
if delta_theta<=theta-theta_m
    flag=3;
    theta=theta-delta_theta;
    x=x2+r2*cos(theta);
    y=y2+r2*sin(theta);
else
    di=sqrt((x-x1)^2+(y-y1)^2);
    delta_theta=acos((di^2+r1^2-d^2)/2/di/r1);
    theta_C1_di=atan((y-y1)/(x-x1));
    theta=theta_C1_di+delta_theta;
    flag=2;
    x=x1+r1*cos(theta);
    y=y1+r1*sin(theta);
end
end

function [x,y,theta,flag]=solve_point_4_3(p,x,y,theta,d,r,xc,yc,theta_m)
a=p/2/pi;
theta=solve_theta2(p,x,y,theta,d);
if theta>=4.5/a-pi
    flag=4;
    x=a*(theta+pi)*cos(theta);
    y=a*(theta+pi)*sin(theta);
else
    fun=@(t)(xc+r*cos(theta_m-t)-x).^2+(yc+r*sin(theta_m-t)-y).^2-d^2;
    q=0.1;
    options = optimoptions('fsolve','Display','off');
    delta_theta=fsolve(fun,theta+q,options);
    theta=theta_m-delta_theta;
    flag=3;
    x=xc+r*cos(theta);
    y=yc+r*sin(theta);
end
end