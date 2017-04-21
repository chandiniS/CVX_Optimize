%Netwons method to minimize the function 
%f(x,y) = exp(x+y) + x^2 + 0.5y^2 +x*y -y st 2*x^2 +3*y^2 <= 1


f_xy = @(x,y)(exp(x+y) + x^2 + 0.5*y^2 + x*y -y ) ; %Original objective function
fb_xy = @(x,y,t)(t*(exp(x+y) + x^2 + 0.5*y^2 + x*y -y) - log(1-2*x^2-3*y^2)); %with log barrier

%Gradient of teh above function
d_fx = @(x,y,t)(t*(exp(x+y)+2*x+y)+ 4*x/(1-2*x^2-3*y^2));
d_fy = @(x,y,t)(t*(exp(x+y)+y+x -1)+ 6*y/(1-2*x^2-3*y^2));
%hessian
d2_fxy = @(x,y,t)([t*(exp(x+y)+2) t*(exp(x+y) +1); t*(exp(x+y) +1) t*(exp(x+y)+1)] + (1/(1-2*x^2-3*y^2)^2) *[4+8*x^2-12*y^2 24*x*y; 24*x*y 6-12*x^2+18*y^2]);
t_val = [1 10 100 1000 10000];
 for tt = t_val
     x_k = 0; y_k=0; xks =[]; yks=[]; fks=[]; iterk =[];
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
       xks = [xks x_k]; yks= [yks y_k]; fks =[fks f_xy(x_k,y_k)];%original objective
       iterk =[iterk iter];
    end
    figure
subplot(3,1,1)
plot(iterk,xks,'LineWidth',1);
xlabel('k');
ylabel('x^{(k)}');

subplot(3,1,2)
plot(iterk,yks,'LineWidth',1);
xlabel('k');
ylabel('y^{(k)}');

subplot(3,1,3)
plot(iterk,fks,'LineWidth',1);
xlabel('k');
ylabel('f(x^{(k)},y^{(k)})');
    
end
