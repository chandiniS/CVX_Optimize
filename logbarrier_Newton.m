%Optimize f: (x-2y)^2/(1+4x+3y) -x -2*y using log barrier method
%Using newton's method for descent step
%Parameter 't' is varied from 1 to 10000

clear
f_xy = @(x,y)((x-2*y)^2 /(1+4*x+3*y)-x-2*y);%Actual function
fb_xy = @(x,y,t)(t*((x-2*y)^2 /(1+4*x+3*y)-x-2*y )-log(1-x^2-y^2)-log(x)-log(y));%with logarithmic barrier
f1 = @(x,y)(1+4*x+3*y);
f2 = @(x,y)(4*x^2-28*y^2+6*x*y+2*x-4*y);
f3 = @(x,y)(-19*x^2+12*y^2+32*x*y-4*x+8*y);
f4 = @(x,y)(1-x^2-y^2);
%Gradient
Dfx =@(x,y,t)(t*(f2(x,y)/f1(x,y)^2 -1) + (2*x)/(f4(x,y)) -1/x);
Dfy =@(x,y,t)(t*(f3(x,y)/f1(x,y)^2 -2) + (2*y)/(f4(x,y)) -1/y);
%Hessian
D2f11 = @(x,y,t)(t*(f1(x,y)*(8*x+6*y+2)-8*f2(x,y))/f1(x,y)^3 +((1+2*x^2-2*y^2)/f4(x,y)^2 )+1/x^2);
D2f12 = @(x,y,t)(t*(f1(x,y)*(-38*x+32*y-4)-8*f3(x,y))/f1(x,y)^3 + (4*x*y)/f4(x,y)^2);
D2f21 = @(x,y,t)(t*(f1(x,y)*(-56*y+6*x-4)-6*f2(x,y))/f1(x,y)^3 + (4*x*y)/f4(x,y)^2);
D2f22 = @(x,y,t)(t*(f1(x,y)*(24*y+32*x+8)-6*f3(x,y))/f1(x,y)^3 +(1+2*y^2-2*x^2)/f4(x,y)^2 +1/y^2);

H= @(x,y,t)([D2f11(x,y,t) D2f12(x,y,t); D2f21(x,y,t) D2f22(x,y,t)]);

ts=[1 10 100 1000 10000];
 x_k = 0.5; y_k=0.5; 
 for tt=10000
     xks =[]; yks=[]; fks=[]; iterk =[];
    alpha = 0.25; beta = 0.6;%intial point
    for iter = 1:100  
        grad = [Dfx(x_k,y_k,tt); Dfy(x_k,y_k,tt)];
        x_nt = -inv(H(x_k,y_k,tt))*grad;%newton step
        lambda = grad' *inv(H(x_k,y_k,tt))*grad;
        if((lambda/2) < 1e-14),  break;     end;
        t =1;
        oldobj = fb_xy(x_k,y_k,tt);
        obj = oldobj;
         while(obj >= oldobj-alpha*t*(x_nt(1)^2 + x_nt(2)^2))
             t = t*beta
             if(f4(x_k+t*x_nt(1),y_k +t*x_nt(2))>0)
                 obj = fb_xy(x_k+t*x_nt(1), y_k +t*x_nt(2),tt)
             end
         end
       x_k = x_k + t*x_nt(1);
       y_k = y_k + t*x_nt(2);   
       xks = [xks x_k]; yks= [yks y_k]; fks =[fks f_xy(x_k,y_k)];%original objective
       iterk =[iterk iter];
    end
    figure
subplot(3,1,1)
plot(iterk,xks,'LineWidth',1,'Marker','*');
xlabel('k');
ylabel('x^{(k)}');

subplot(3,1,2)
plot(iterk,yks,'LineWidth',1,'Marker','*');
xlabel('k');
ylabel('y^{(k)}');

subplot(3,1,3)
plot(iterk,fks,'LineWidth',1,'Marker','*');
xlabel('k');
ylabel('f(x^{(k)},y^{(k)})');
 end

