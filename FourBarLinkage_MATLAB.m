%FOUR BAR LINKAGE MECHANISM
%BY HUSSEIN SHATA

clear all
close all

L0 = 7;
L1 = 3;
L2 = 7.5;
L3 = 4.5;

L_AP = 4;
alpha = 35;

theta2 = 0:2:360;
figure(1); 
for i = 1:length(theta2)
    
  if theta2(i) > 360
   theta2(i) = theta2(i)-360;
  end
  
   AC(i) = sqrt(L0^2 + L1^2 - 2*L0*L1*cosd(theta2(i)));
   beta(i) = acosd((L0^2 + AC(i)^2 - L1^2) / (2*L0*AC(i)));
   psi(i) = acosd((L2^2 + AC(i)^2 - L3^2) / (2*L2*AC(i)));
   lamda(i) = acosd((L3^2 + AC(i)^2 - L2^2) / (2*L3*AC(i)));
   
   theta3(i) = psi(i) - beta(i);
   theta4(i) = 180 - lamda(i) - beta(i);
   
   if theta2(i) > 180
       theta3(i) = psi(i) + beta(i);
       theta4(i) = 180 - lamda(i) + beta(i);
   end
   
   %Define the joints position
   Ox(i) = 0;
   Oy(i) = 0;
   
   Ax(i) = Ox(i) + L1*cosd(theta2(i));
   Ay(i) = Oy(i) + L1*sind(theta2(i));
   
   Bx(i) = Ox(i) + Ax(i) + L2*cosd(theta3(i));
   By(i) = Oy(i) + Ay(i) + L2*sind(theta3(i));
   
   Cx(i) = L0;
   Cy(i) = 0;
   
   Px(i) = Ox(i) + Ax(i) + L_AP*cosd(alpha+theta3(i));
   Py(i) = Oy(i) + Ay(i) + L_AP*sind(alpha+theta3(i));
   
   
   r = [ L2*cosd(theta3(i)), -L3*cosd(theta4(i));...
         L2*sind(theta3(i)), -L3*sind(theta4(i)) ];
    
   v = [ -L1*alpha*cosd(theta2(i)); -L1*alpha*sind(theta2(i))];
   
   w(i,:) = inv(r)*v;
   om2(i) = w(i,1);
   om4(i) = w(i,2);
   
   alpha_dot(i) = - (L1*om2(i)*cosd(theta2(i)) - L3*om4(i)*cosd(theta4(i)))/ ...
                    (L2*cosd(alpha) + L_AP*cosd(alpha));
                
   V_Px(i) = alpha_dot(i)*Px(i);
   V_Py(i) = alpha_dot(i)*Py(i);

   plot( [Ox(i) Ax(i)], [Oy(i) Ay(i)], [Ax(i) Bx(i)], [Ay(i) By(i)] ...
       , [Bx(i) Cx(i)], [By(i) Cy(i)], 'LineWidth',3)
   hold on;
   
   plot( [Ax(i) Px(i)], [Ay(i) Py(i)],'LineWidth',3)
   plot( [Bx(i) Px(i)], [By(i) Py(i)],'LineWidth',3)
  
    grid;
    axis equal;
    axis([-5 15 -5 10]);
    drawnow;
    hold off;
    
end

figure(2)
plot(theta2, om2,'LineWidth',2)
title('Angular Velocity of \theta_3')
xlabel('\theta_2')
ylabel('Joint AB')

figure(3)
plot(theta2, om4,'LineWidth',2)
title('Angular Velocity of \theta_4')
xlabel('\theta_2')
ylabel('Joint AB')

figure(4)
plot(theta2,  V_Px,'LineWidth',2)
title('x-component velocity of P')
xlabel('\theta_2')
ylabel('Velocity of Point P')

figure(5)
plot(theta2,  V_Py,'LineWidth',2)
title('y-component velocity of P')
xlabel('\theta_2')
ylabel('Velocity of Point P')