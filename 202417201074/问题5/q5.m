p=1.7; % 螺距
a = p / (2 * pi); % 螺线系数
d1 = 2.86; % 龙头两把手间距
d2 = 1.65; % 除龙头板外其他板两把手间距
V_0=1:0.1:2; % 龙头运动速度
V_max=zeros(1,length(V_0));

for v=1:length(V_0)
    v0=V_0(v);
    R=4.5;
    theta_ru=R/a;
    theta_chu=R/a-pi;
    slope=(a*sin(theta_ru)+R*cos(theta_ru))/(a*cos(theta_ru)-R*sin(theta_ru));
    theta_max_a=atan(-1/slope)+pi;
    theta_d=atan(tan(theta_ru))+pi-theta_max_a;
    rC1_C2=R/cos(theta_d);
    rc2=rC1_C2/3;
    rc1=rc2*2;
    phi=2*theta_d;
    Sc1=rc1*(pi-phi);
    Sc2=rc2*(pi-phi);
    theta_min_a=theta_max_a-Sc1/rc1;
    theta_min_b=theta_min_a-pi;
    theta_max_b=theta_min_b+Sc2/rc2;
    % 得到两个圆圆心坐标
    x_C1=R*cos(theta_ru)+rc1*cos(theta_max_a-pi);
    y_C1=R*sin(theta_ru)+rc1*sin(theta_max_a-pi);
    x_C2=R*cos(theta_chu)-rc2*cos(theta_max_b);
    y_C2=R*sin(theta_chu)-rc2*sin(theta_max_b);
    dt=0.1;
    Ttotal=(Sc1+Sc2)/v0;
    T=0:dt:Ttotal;
    X=nan*zeros(224,length(T));
    Y=nan*zeros(224,length(T));
    T=nan*zeros(224,length(T));
    tt_c1=0:dt:Sc1/v0;
    theta_C1=-v0*tt_c1/rc1+theta_max_a;
    T(1,(1:length(tt_c1)))=theta_C1;
    X(1,(1:length(tt_c1)))=rc1*cos(theta_C1)+x_C1;
    Y(1,(1:length(tt_c1)))=rc1*sin(theta_C1)+y_C1;
    tt_c2=tt_c1(end)+dt:dt:(Sc1+Sc2)/v0;
    theta_C2=v0*(tt_c2-Sc1/v0)/rc2+theta_min_a-pi;
    T(1,length(theta_C1)+(1:length(tt_c2)))=theta_C2;
    X(1,length(theta_C1)+(1:length(tt_c2)))=rc2*cos(theta_C2)+x_C2;
    Y(1,length(theta_C1)+(1:length(tt_c2)))=rc2*sin(theta_C2)+y_C2;
    
    % 循环求解其余板凳位置
    N=223; % 板凳总个数
    ttotal=0:dt:(Sc1+Sc2)/v0;
    for jj=(1:length(ttotal))
        j=jj;
        if ttotal(jj)<0 
            for i=2:N+1
                d=d1*(i<=2)+d2*(i>2);
                thetaij=solve_theta1(p,X(i-1,j),Y(i-1,j),T(i-1,j),d);
                T(i,j)=thetaij;
                X(i,j)=a*thetaij*cos(thetaij);
                Y(i,j)=a*thetaij*sin(thetaij);
            end
        elseif ttotal(jj)>=0 && ttotal(jj)<=Sc1/v0 
            flag=2;
            for i=2:N+1
                d=d1*(i<=2)+d2*(i>2);
                if flag==2
                    [xi,yi,thetai,flag]=solve_point_2_1(p,X(i-1,j),Y(i-1,j),T(i-1,j),d,rc1,x_C1,y_C1,theta_max_a);
                    T(i,j)=thetai;
                    X(i,j)=xi;
                    Y(i,j)=yi;
                else
                    thetaij=solve_theta1(p,X(i-1,j),Y(i-1,j),T(i-1,j),d);
                    T(i,j)=thetaij;
                    X(i,j)=a*thetaij*cos(thetaij);
                    Y(i,j)=a*thetaij*sin(thetaij);
                end
            end
        elseif ttotal(jj)>Sc1/v0 && ttotal(jj)<=Sc1/v0+Sc2/v0
            flag=3;
            for i=2:N+1
                d=d1*(i<=2)+d2*(i>2);
                if flag==3
                    [xi,yi,thetai,flag]=solve_point_3_2(X(i-1,j),Y(i-1,j),T(i-1,j),d,rc1,x_C1,y_C1,rc2,x_C2,y_C2,theta_min_b);
                    T(i,j)=thetai;
                    X(i,j)=xi;
                    Y(i,j)=yi;
                elseif flag==2
                    [xi,yi,thetai,flag]=solve_point_2_1(p,X(i-1,j),Y(i-1,j),T(i-1,j),d,rc1,x_C1,y_C1,theta_max_a);
                    T(i,j)=thetai;
                    X(i,j)=xi;
                    Y(i,j)=yi;
                else
                    thetaij=solve_theta1(p,X(i-1,j),Y(i-1,j),T(i-1,j),d);
                    T(i,j)=thetaij;
                    X(i,j)=a*thetaij*cos(thetaij);
                    Y(i,j)=a*thetaij*sin(thetaij);
                end
            end
        else
            flag=4;
            for i=2:N+1
                d=d1*(i<=2)+d2*(i>2);
                if flag==4 
                    [xi,yi,thetai,flag]=solve_point_4_3(p,X(i-1,j),Y(i-1,j),T(i-1,j),d,rc2,x_C2,y_C2,theta_max_b);
                    T(i,j)=thetai;
                    X(i,j)=xi;
                    Y(i,j)=yi;
                elseif flag==3
                    [xi,yi,thetai,flag]=solve_point_3_2(X(i-1,j),Y(i-1,j),T(i-1,j),d,rc1,x_C1,y_C1,rc2,x_C2,y_C2,theta_min_b);
                    T(i,j)=thetai;
                    X(i,j)=xi;
                    Y(i,j)=yi;
                elseif flag==2 
                    [xi,yi,thetai,flag]=solve_point_2_1(p,X(i-1,j),Y(i-1,j),T(i-1,j),d,rc1,x_C1,y_C1,theta_max_a);
                    T(i,j)=thetai;
                    X(i,j)=xi;
                    Y(i,j)=yi;
                else 
                    thetaij=solve_theta1(p,X(i-1,j),Y(i-1,j),T(i-1,j),d);
                    T(i,j)=thetaij;
                    X(i,j)=a*thetaij*cos(thetaij);
                    Y(i,j)=a*thetaij*sin(thetaij);
                end
            end
        end
        waitbar((length(ttotal)*(v-1)+jj)/(length(ttotal)*length(V_0)),hwait,'已经完成...')
    end
    
    Vx=zeros(size(X));
    Vy=Vx;
    Vx(:,1)=(X(:,2)-X(:,1))/dt;
    Vx(:,end)=(X(:,end)-X(:,end-1))/dt;
    Vx(:,2:end-1)=(X(:,3:end)-X(:,1:end-2))/2/dt;
    Vy(:,1)=(Y(:,2)-Y(:,1))/dt;
    Vy(:,end)=(Y(:,end)-Y(:,end-1))/dt;
    Vy(:,2:end-1)=(Y(:,3:end)-Y(:,1:end-2))/2/dt;
    V=sqrt(Vx.^2+Vy.^2);
    V_max(v)=max(max(V));
    
end

% 下面可视化
p=polyfit(V_0,V_max,1); % 线性拟合得到Vmax与v0的解析关系
figure
p1=plot(V_0,V_max,'ro','LineWidth',1);
xlabel('龙头速度')
ylabel('队伍把手最高速度')
hold on
color=rand(1,3);
p2=plot(linspace(V_0(1),V_0(end),100),polyval(p,linspace(V_0(1),V_0(end),100)),'Color',color);
plot([V_0(1) V_0(end)],[2 2],'k--')
legend('队伍最大速度与龙头速度数据点','线性拟合关系')
v0max=(2-p(2))/p(1);
disp(['龙头最大行进速度为:',num2str(v0max)])

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

function [x,y,theta,flag]=solve_point_2_1(p,x1,y1,theta,d,rC1,x,y,theta_m)
a=p/2/pi;
delta_theta=2*asin(d/2/rC1);
if delta_theta<=theta_m-theta
    flag=2;
    theta=theta+delta_theta;
    x=x+rC1*cos(theta);
    y=y+rC1*sin(theta);
else
    theta=solve_theta1(p,x1,y1,4.5/a,d);
    flag=1;
    x=a*theta*cos(theta);
    y=a*theta*sin(theta);
end
end

function [x,y,theta,flag]=solve_point_3_2(x,y1,theta1,d,rC1,x_c1,y_c1,rC2,x2,y2,theta_m)
delta_theta=2*asin(d/2/rC2);
if delta_theta<=theta1-theta_m
    flag=3;
    theta=theta1-delta_theta;
    x=x2+rC2*cos(theta);
    y=y2+rC2*sin(theta);
else
    di=sqrt((x-x_c1)^2+(y1-y_c1)^2);
    delta_theta=acos((di^2+rC1^2-d^2)/2/di/rC1);
    theta_C1_di=atan((y1-y_c1)/(x-x_c1));
    theta=theta_C1_di+delta_theta;
    flag=2;
    x=x_c1+rC1*cos(theta);
    y=y_c1+rC1*sin(theta);
end
end

function [x,y,theta,flag]=solve_point_4_3(p,x,y,theta,d,rc2,xc,yc,theta_m)
a=p/2/pi;
theta=solve_theta2(p,x,y,theta,d);
if theta>=4.5/a-pi 
    x=a*(theta+pi)*cos(theta);
    y=a*(theta+pi)*sin(theta);
else
    fun=@(t)(xc+rc2*cos(theta_m-t)-x).^2+(yc+rc2*sin(theta_m-t)-y).^2-d^2;
    q=-0.1;
    options = optimoptions('fsolve','Display','off');
    delta_theta=fsolve(fun,theta_m+q,options);
    theta=theta_m-delta_theta;
    flag=3;
    x=xc+rc2*cos(theta);
    y=yc+rc2*sin(theta);
end
end