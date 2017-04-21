%Optimize e^(x+y) + x2+ 0.5y2 +xy -y with Gradient descent & ELS

clear

f_xy = @(x,y)(exp(x+y)+x^2+0.5*y^2+ x*y -y);
df_x = @(x,y)(exp(x+y)+2*x+y);
df_y = @(x,y)(exp(x+y)+y+x-1);

x_k =1; y_k = 1;
    t_els = 0; p_els = inf;
    %negative gradient
    Dx = -df_x(x_k,y_k); Dy = -df_y(x_k,y_k);
    
    
obj_t =[];t_x=0:0.001:0.99;
    for t = 0:0.001:0.99
         obj = f_xy(x_k+t*(Dx),y_k+t*(Dy))% f(x+t*dx)
         obj_t= [obj_t obj];
         if (obj<p_els)
             t_els =t
             p_els=obj;
         end        
    end
figure
plot(t_x,obj_t,'LineWidth',1);
xlabel('t^{(0)}');
ylabel('Objective value- f(x^{(0)}+t^{(0)}\Deltax^{(0)}, y^{(0)}+t^{(0)}\Deltay^{(0)})');