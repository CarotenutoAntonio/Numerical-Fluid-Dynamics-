clc;
clear;
close all;


I_K0Re2000=load("I_K0Re2000_vero.mat");

I_K03Re2000=load("I_K03Re2000.mat");

I_K2Re2000=load("I_K2Re2000.mat");

I_K5Re2000=load("I_K5Re2000.mat");

I_K1Re2000=load("I_K1Re2000_vero.mat");

n=round((3/5)*length(I_K1Re2000.I));
I_K1Re2000.I=I_K1Re2000.I(1:n);
I_K1Re2000.time=I_K1Re2000.time(1:n);
I_K1Re500=load("I_K1Re500_vero.mat");

I_K1Re1000=load("I_K1Re1000_vero.mat");

I_K1Re5000=load("I_K1Re5000_vero.mat");

I_omega05Re2000=load("I_omega05Re5000_vero.mat");

I_omega08Re2000=load("I_omega08Re2000_vero.mat");

I_omega12Re2000=load("I_omega12Re2000_vero.mat");

I_omega2Re2000=load("I_omega2Re5000_vero.mat");


   figure(1);
   plot(I_K0Re2000.time,I_K0Re2000.I, 'r',I_K03Re2000.time,I_K03Re2000.I, 'g',...
       I_K1Re2000.time,I_K1Re2000.I, 'y',I_K2Re2000.time,I_K2Re2000.I, 'b',...
       I_K5Re2000.time,I_K5Re2000.I, 'k'); 
   hold on;
   legend('k=0','k=0.3','k=1','k=2','k=5','Location','northwest')
   title('Particle dispersion, effect of k');
   xlabel('t'); ylabel('I/I_0'); drawnow;



    figure(2);
   plot(I_K1Re500.time,I_K1Re500.I, 'r',I_K1Re1000.time,I_K1Re1000.I, 'g',...
       I_K1Re2000.time,I_K1Re2000.I, 'y',I_K1Re5000.time,I_K1Re5000.I,  'k'); 
   hold on;
   legend('Re=500','Re=1000','Re=2000','Re=5000','Location','northwest')
   title('Particle dispersion, effect of Re');
   xlabel('t'); ylabel('I/I_0'); drawnow;


       figure(3);
   plot(I_omega05Re2000.time,I_omega05Re2000.I, 'r',I_omega08Re2000.time,I_omega08Re2000.I, 'g',...
       I_K1Re2000.time,I_K1Re2000.I,'y',I_omega12Re2000.time,I_omega12Re2000.I, 'b',...
       I_omega2Re2000.time,I_omega2Re2000.I,  'k'); 
   hold on;
   legend('\omega_{ref} =0.5','\omega_{ref} =0.8','\omega_{ref} =1','\omega_{ref} =1.2',...
       '\omega_{ref} =2','Location','northwest');

   title('Particle dispersion, effect of $\omega_{ref}$ ','Interpreter','latex');
   xlabel('t'); ylabel('I/I_0'); drawnow;

