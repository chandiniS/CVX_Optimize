%Optimize e^(x+y) + x2+ 0.5y2 +xy -y with Gradient descent & BLS
clear

f_xy = @(x,y)(exp(x+y)+x^2+0.5*y^2+ x*y -y);
df_x = @(x,y)(exp(x+y)+2*x+y);
df_y = @(x,y)(exp(x+y)+y+x-1);

x_k =1; y_k = 1;   
 obj_k =[];t_k=[];  xks=[]; yks=[]; 
for iter = 1:50 
     t_els = 0; p_els = inf;
      %negative gradient
    Dx = -df_x(x_k,y_k); Dy = -df_y(x_k,y_k);
      
 if(norm([Dx Dy]) < 1e-13)
    break;
end
    %els
    for t = 0:0.01:10
         obj = f_xy(x_k+t*(Dx),y_k+t*(Dy))% f(x+t*dx)
        
         if (obj<p_els)
             t_els =t
             p_els=obj;
         end        
    end
    
    x_k = x_k+t_els*Dx;
    y_k = y_k+t_els*Dy;
    xks =[xks x_k];
    yks=[yks y_k];
    obj_k= [obj_k f_xy(x_k,y_k)];
    t_k = [t_k t_els];
end
iter= 1:1:50;
figure
subplot(4,1,1)
plot(iter,xks,'LineWidth',1);
xlabel('k');
ylabel('x^{(k)}');

subplot(4,1,2)
plot(iter,yks,'LineWidth',1);
xlabel('k');
ylabel('y^{(k)}');

subplot(4,1,3)
plot(iter,t_k,'LineWidth',1);
xlabel('k');
ylabel('t^{(k)}');

subplot(4,1,4)
plot(iter,obj_k,'LineWidth',1);
xlabel('k');
ylabel('f(x^{(k)},y^{(k)})');