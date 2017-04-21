%Optimize e^(x+y) + x2+ 0.5y2 +xy -y with Newton's Method & BLS
clear
f_xy = @(x,y)(exp(x+y)+x^2+0.5*y^2+ x*y -y);
df_x = @(x,y)(exp(x+y)+2*x+y);%gradient
df_y = @(x,y)(exp(x+y)+y+x-1);%gradient
d2f = @(x,y)([exp(x+y)+2 exp(x+y)+1; exp(x+y)+1 exp(x+y)+1]);%Hessian
x_k = 1; y_k=1; obj_k = []; t_k =[]; xks= []; yks=[];iterk=[];
for iter = 1:100
    xnt = -inv(d2f(x_k,y_k)) * [df_x(x_k,y_k); df_y(x_k,y_k)];
    Dx =xnt(1); Dy = xnt(2);  
    lambda2= [df_x(x_k,y_k) df_y(x_k,y_k)]* inv(d2f(x_k,y_k)) * [df_x(x_k,y_k); df_y(x_k,y_k)];   
    if(lambda2/2 <= 1e-9),break;   end;   
    t = 1; alpha = 1/2; beta = 1/2;    
    oldobj = f_xy(x_k,y_k);
    obj=oldobj;
         while(obj>=(oldobj - alpha*t*(Dx^2+Dy^2)))
            t = t* beta;
            obj = f_xy(x_k+t*Dx, y_k+t*Dx)          
         end       
    x_k = x_k+t*Dx;
    y_k = y_k+t*Dy;
    
    xks =[xks x_k];    yks=[yks y_k];   obj_k= [obj_k f_xy(x_k,y_k)];  t_k = [t_k t]; iterk= [iterk iter];   
end
figure
subplot(4,1,1)
plot(iterk,xks,'LineWidth',1);
xlabel('k');
ylabel('x^{(k)}');

subplot(4,1,2)
plot(iterk,yks,'LineWidth',1);
xlabel('k');
ylabel('y^{(k)}');

subplot(4,1,3)
plot(iterk,t_k,'LineWidth',1);
xlabel('k');
ylabel('t^{(k)}');

subplot(4,1,4)
plot(iterk,obj_k,'LineWidth',1);
xlabel('k');
ylabel('f(x^{(k)},y^{(k)})');