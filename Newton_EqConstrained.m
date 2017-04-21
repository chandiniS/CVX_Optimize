%Sove e^(x+y) +x2 +0.5y2 +xy -y such that x+2y =3
%Eliminate equality constraints  and reformulate as unconstrained
%minimization- Newton's Method
clear
f_xy = @(x,y)(exp(x+y)+x^2+ 0.5*y^2+x*y -y);
A = [1 2];
df_x = @(x,y)(exp(x+y)+2*x+y);
df_y = @(x,y)(exp(x+y)+y+x-1);
d2f = @(x,y)([exp(x+y)+2 exp(x+y)+1; exp(x+y)+1 exp(x+y)+1]);
x_k = 1; y_k=1; obj_k = []; t_k =[]; xks= []; yks=[];iterk=[];
for iter = 1:100
    xnt = inv([d2f(x_k,y_k) A'; A 0]) * [ -df_x(x_k,y_k); -df_y(x_k,y_k); 0];
    vx = xnt(1); vy = xnt(2);   
    lambda_x2 = [vx vy] * d2f(x_k,y_k) * [vx; vy];    
    if(lambda_x2/2< 1e-12),       break; end;
    t = 1; alpha = 1/4; beta = 1/2;   
    oldobj = f_xy(x_k,y_k);
    obj=oldobj
         while(obj>=(oldobj - alpha*t*([df_x(x_k,y_k) df_y(x_k,y_k)]*[vx;vy])))
            t = t* beta
            obj = f_xy(x_k+t*vx, y_k+t*vy)          
         end        
    x_k = x_k+t*vx;
    y_k = y_k+t*vy;
    
     xks =[xks x_k]; yks=[yks y_k]; obj_k= [obj_k f_xy(x_k,y_k)]; t_k = [t_k t]; iterk= [iterk iter];
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