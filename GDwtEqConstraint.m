%Sove e^(x+y) +x2 +0.5y2 +xy -y such that x+2y =3
%Eliminate equality constraints  and reformulate as unconstrained
%minimization
clear

f_y = @(y)(exp(3-y)+2.5*y^2 -10*y +9);
df_y = @(y)(-exp(3-y)+5*y-10);
x_k =1; y_k = 1;   
 obj_k =[];t_k=[];  xks=[]; yks=[]; iterk=[];
for iter = 1:100 
    Dy = -df_y(y_k);
     if(norm([Dy]) < 1e-5)  break;   end
    %BLS
    alpha =1/12; beta = 0.9; t=1; 
         oldobj = f_y(y_k);
         obj=oldobj;
         while(obj>=(oldobj - alpha*t*Dy^2))
            t = t* beta
            obj = f_y(y_k+t*Dy)           
         end       
    x_k = 3 - 2*y_k;
    y_k = y_k+t*Dy
    xks =[xks x_k]; yks=[yks y_k]; obj_k= [obj_k f_y(y_k)]; t_k = [t_k t]; iterk= [iterk iter];
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