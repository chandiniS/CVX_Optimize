%%Newton's method to minimize the function. Obtain the central path
%f(x,y) = exp(x+y) + x^2 + 0.5y^2 +x*y -y st 2*x^2 +3*y^2 <= 1
f_xy = @(x,y)(exp(x+y) + x^2 + 0.5*y^2 + x*y -y ) ;
fb_xy = @(x,y,t)(t*(exp(x+y) + x^2 + 0.5*y^2 + x*y -y) - log(1-2*x^2-3*y^2)); %with barrier

%Gradient of teh above function
d_fx = @(x,y,t)(t*(exp(x+y)+2*x+y)+ 4*x/(1-2*x^2-3*y^2));
d_fy = @(x,y,t)(t*(exp(x+y)+y+x -1)+ 6*y/(1-2*x^2-3*y^2));
%hessian
d2_fxy = @(x,y,t)([t*(exp(x+y)+2) t*(exp(x+y) +1); t*(exp(x+y) +1) t*(exp(x+y)+1)] + (1/(1-2*x^2-3*y^2)^2) *[4+8*x^2-12*y^2 24*x*y; 24*x*y 6-12*x^2+18*y^2]);
t_val = [1 10 100 1000];
xts =[]; yts=[]; 
 for tt = 1:10000
     x_k = 0; y_k=0; 
    alpha = 0.5; beta = 0.9;%intial point
    for iter = 1:100  
        grad = [d_fx(x_k,y_k,tt);d_fy(x_k,y_k,tt)];
        x_nt = -inv(d2_fxy(x_k,y_k,tt))*grad;%newton step
        lambda = grad' *inv(d2_fxy(x_k,y_k,tt)) *grad;
        if((lambda/2) < 1e-12),  break;     end;
        t =1;
        oldobj = fb_xy(x_k,y_k,tt);
        obj = oldobj;
         while(obj >= oldobj-alpha*t*(x_nt(1)^2 + x_nt(2)^2))
             t = t*beta
             if(x_k+t*x_nt(1) + y_k +t*x_nt(2) <=1)
                 obj = fb_xy(x_k+t*x_nt(1), y_k +t*x_nt(2),tt);
             end
         end
       x_k = x_k + t*x_nt(1);
       y_k = y_k + t*x_nt(2);          
    end
   xts=[xts x_k];
   yts=[yts y_k];
 end
 
figure
plot(xts,yts,'LineWidth',1,'Marker','*');
xlabel('x^{*}(t)');
ylabel('y^{*}(t)');
    

