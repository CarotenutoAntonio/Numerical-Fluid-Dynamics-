clear; close all; clc;
% Esercitazione dal corso di Fluidodinamica Numerica
% Prof. G. Coppola
%  
% Codice di esempio per la discretizzazione delle equazioni di
% Navier-Stokes 2D in variabili primitive u,v,p.
% La procedura di integrazione segue il classico metodo di proiezione.
% La discretizzazione spaziale si basa sullo staggering alla Harlow-Welch, 
% nel quale le variabili componenti di velocità e pressione sono 
% localizzate sul mesh secondo lo schema:
% 
% 
%                                V_(i,j+1)
%                            o-----^-----o
%                            |           |
%                            |   p_(i,j) |
%                    U_(i,j) >     x     > U_(i+1,j)
%                            |           |
%                            |           |
%                            o-----^-----o
%                                V_(i,j)
% 
%  con i simboli:
%     x centro-cella (pressione)
%     o vertici di cella
%     > facce verticali (componenti u)
%     ^ facce orizzontali (componenti v)
%
% Il numero di incognite per U e V è dato da:
%
%                   U  --> Nx*(Ny-1)     
%                   V  --> Ny*(Nx-1)     
%                   p  --> (Nx-1)*(Ny-1) 

% 
% La integrazione temporale viene effettuata mediate il classico schema 
% RK4, nel quale la proiezione del campo di velocità sul sottospazio a
% divergenza nulla viene fatta ad ogni stage.
% 

global Nx Ny Re h hq LapPsi hx hy x y Lx Ly Udwn Uup         % Dichiarazione globale delle varabili
global UNord USud UWest UEst VNord VSud VWest VEst
warning off
Lx = 1;    Ly = 1;                   
% Dimensioni del dominio

Nx = 80;   Ny = 80;   N = Nx;        % Numero di nodi lungo ogni lato
x  = linspace(0,Lx,Nx);              % Mesh (uniforme) lungo x
y  = linspace(0,Ly,Ny);              % Mesh (uniforme) lungo y



hx = x(2) - x(1);                    % Passo spaziale lungo x
hy = y(2) - y(1);                    % Passo spaziale lungo y
h  = hx;           hq = h*h;        
Re = 2000;                           % Numero di Reynolds
T  = 30;  


%% modifico griglia per u e v
x_int_U = x;                         % Mesh per interpolazione nei nodi U
y_int_U = y(1:end-1)+hy/2;    

x_int_V = x(1:end-1)+hx/2;           % Mesh per interpolazione nei nodi V
y_int_V = y;

[y_meshU, x_meshU] = meshgrid(y_int_U ,x_int_U);
[y_meshV, x_meshV] = meshgrid(y_int_V ,x_int_V);
%%
% Tempo finale della simulazione
% Preallocazione delle variabili. 
% Le variabili U e V sono preallocate in array che contengono anche le
% variabili di bordo che cadono sul boundary. Per la variabile P si
% considerano invece solo i valori interni al dominio.
% La istruzione di preallocazione introduce anche le condizioni iniziali.

U  = zeros(Nx,Ny-1);     V = zeros(Nx-1,Ny);    P = zeros(Nx-1,Ny-1);


Uref  = 1;                           % Velocità di riferimento
Omega_ref=0.8;                         %pulsazione di riferimento
Om=Omega_ref;                        %pulasazione della parete superiore 
K=1;                               %rapporto tra pulasazione inferiore e superiore

C     = 1.2;       beta = 0.8;       % Parametri della discretizzazione numerica
Dt    = min([C*h/Uref,beta*hq*Re]);  % Dt per la stabiità
Nt    = round(T/Dt);                 % Numero di step temporali
time=1:Nt;    time=time*Dt;
I=zeros(1,Nt);                     % Area occupata dalle particelle              


%definiamo la posizione iniziale delle paticelle di polvere depositate sul
%fondo del contenitore 

x_0=[  Lx/2         Ly/2
       Lx/2         9*Ly/20
       Lx/2         11*Ly/20
       9*Lx/20      Ly/2
       9*Lx/20      9*Ly/20
       9*Lx/20      11*Ly/20
       11*Lx/20     Ly/2
       11*Lx/20     9*Ly/20
       11*Lx/20     11*Ly/20];

ColorDots=[1   0   0
           0   1   0
           0   0   1
           0   1   1
           1   0   1
           1   1   0
           0   0   0
           0.4940  0.1840  0.5560
           0.929   0.6940  0.1250];

figure(7);

s1=scatter(x_0(:,1), x_0(:,2), 20, ColorDots,'filled' );
axis([0 Lx 0 Ly]); axis square;

x_d=x_0;
Centroid=[ sum(x_d(:,1)), sum(x_d(:,2))];  Centroid=Centroid/length(x_d(:,1));
Dr=[ x_d(:,1)-Centroid(1)  , x_d(:,2)-Centroid(2)]; Dr2=Dr(:,1).^2+Dr(:,2).^2;
I0=sum(Dr2)/length(x_d(:,1));
   

%imponiamo le condizioni al contorno sulle pareti
UWest = 0; UEst = 0; VNord = 0; VSud = 0;
UNord = Uref*sin(Om*time);  USud = Uref*cos(K*Om*time) ; VWest = 0; VEst = 0;

% Inserimento delle BCs negli array U e V. Questi valori di bordo non
% saranno aggiornati durante la integrazione.
U(1,:)  = UWest;     U(Nx,:) = UEst; 
V(:,Ny) = VNord;     V(:,1)  = VSud;
% Calcolo dell'operatore di Laplace per l'ellittica di pressione
G = numgrid('S',N+1);    Lap = -delsq(G)/hq;
% Implementazione delle condizioni al contorno alla Neumann sull'operatore
% di Laplace discreto
for i = 2:Nx
    j = 2;    k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
    j = Ny;   k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
end
for j = 2:Ny
    i = 2;    k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
    i = Nx;   k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
end
% Calcolo dell'operatore di Laplace con condizioni al contorno alla 
% Dirichlet su un mesh collocato, con nodi ai vertici del nostro mesh.
% Questo operatore verrà usato per calcolare la Psi da U e V ai fini della
% grafica.
G   = numgrid('S',Nx);     LapPsi = -delsq(G);     LapPsi = LapPsi/hq;
% Avanzamento temporale
for it = 1:Nt

    Uup=UNord(it); Udwn=USud(it);
%%%%%%%%%%%%%%%%
% STAGE 1 
    U1 = U;         V1 = V;   xd1=x_d;         [Fu1,Fv1] = RHS_HW_OL(U1,V1,Uup, Udwn);

    %%modifico l'interpolazione
    U_xd1 = interp2(y_meshU, x_meshU, U1, xd1(1:end,2),xd1(1:end,1), 'spline');
    V_xd1 = interp2(y_meshV, x_meshV, V1, xd1(1:end,2),xd1(1:end,1), 'spline');
    UV_int1=[U_xd1(:),V_xd1(:)];
%%%%%%%%%%%%%%%%
% STAGE 2
% Calcolo del campo asteriscato
    U2  = U + Dt*0.5*Fu1;             V2 = V + Dt*0.5*Fv1;
% Soluzione della ellittica di pressione
    DIV = DivCalc(U2,V2);             P = Lap\DIV(:);  
% Calcolo del gradiente di pressione    
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py]   = GradCalc(P);
% Step di proiezione
    U2  = U2 - Px;   V2 = V2 - Py;   [Fu2,Fv2] = RHS_HW_OL(U2,V2,Uup, Udwn);
% Calcolo della posizione delle particelle
    xd2= x_d + Dt*0.5*UV_int1;
%Eseguiamo una nuova interpolazione di U e V sui punti  di xd2

    U_xd2 = interp2(y_meshU, x_meshU, U2, xd2(1:end,2),xd2(1:end,1), 'spline');
    V_xd2 = interp2(y_meshV, x_meshV, V2, xd2(1:end,2),xd2(1:end,1), 'spline');
    UV_int2=[U_xd2(:),V_xd2(:)];

%%%%%%%%%%%%%%%%    
% STAGE 3
    U3  = U + Dt*0.5*Fu2;             V3 = V + Dt*0.5*Fv2;
    DIV = DivCalc(U3,V3);             P  = Lap\DIV(:);   
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py]   = GradCalc(P);
    U3  = U3 - Px;   V3 = V3 - Py;   [Fu3,Fv3] = RHS_HW_OL(U3,V3,Uup, Udwn);

    % Calcolo della posizione delle particelle
    xd3= x_d + Dt*0.5*UV_int2;
    %Eseguiamo una nuova interpolazione di U e V sui punti  di xd3
    U_xd3 = interp2(y_meshU, x_meshU, U3, xd3(1:end,2),xd3(1:end,1), 'spline');
    V_xd3 = interp2(y_meshV, x_meshV, V3, xd3(1:end,2),xd3(1:end,1), 'spline');
    UV_int3=[U_xd3(:),V_xd3(:)];

%%%%%%%%%%%%%%%%    
% STAGE 4
    U4  = U + Dt*Fu3;                 V4 = V + Dt*Fv3;
    DIV = DivCalc(U4,V4);             P  = Lap\DIV(:);   
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py] = GradCalc(P);
    U4  = U4 - Px;   V4 = V4 - Py;   [Fu4,Fv4] = RHS_HW_OL(U4,V4,Uup, Udwn);

    % Calcolo della posizione delle particelle
    xd4= x_d + Dt*UV_int3;
%Eseguiamo una nuova interpolazione di U e V sui punti  di xd2
    U_xd4 = interp2(y_meshU, x_meshU, U4, xd4(1:end,2),xd4(1:end,1), 'spline');
    V_xd4 = interp2(y_meshV, x_meshV, V4, xd4(1:end,2),xd4(1:end,1), 'spline');
    UV_int4=[U_xd4(:),V_xd4(:)];

%%%%%%%%%%%%%%%%    
%   Step Finale
    U   = U + Dt*((1/6)*(Fu1 + Fu4) + (1/3)*(Fu2 + Fu3));    
    V   = V + Dt*((1/6)*(Fv1 + Fv4) + (1/3)*(Fv2 + Fv3));
    DIV = DivCalc(U,V);               P  = Lap\DIV(:);
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py] = GradCalc(P);
    U   = U - Px;    V = V - Py;

    % Calcolo della posizione delle particelle
    x_d = x_d + Dt*((1/6)*(UV_int1 + UV_int4) + (1/3)*(UV_int2 + UV_int3));
    
    %The following routine checks that no particle excedes the boundaries
    %of the domain. 

    for i=1:9
        if x_d(i,1)<0
            x_d(i,1)=0;
        elseif x_d(i,1)>Lx
            x_d(i,1)=Lx;
        end
        if x_d(i,2)<0
            x_d(i,2)=0;
        elseif x_d(i,2)>Ly
            x_d(i,2)=Ly;
        end
    end


   %%%%%%%%%%%%%%%%%%%%%%%%

   %Let's draw the Pathlines
   
   figure(1);
   
   s1=scatter(x_d(:,1), x_d(:,2), 20, ColorDots ,'filled');
   hold on
   axis([0 Lx 0 Ly]); axis square;
   title(['Particle scattering t= ', num2str(it*Dt)]); drawnow;
  
   %Consideriamo la dispersione delle paicelle rispetto al centroide
   %definito come segue 
   Centroid=[ sum(x_d(:,1)), sum(x_d(:,2))];  Centroid=Centroid/length(x_d(:,1));
   
   Dr=[ x_d(:,1)-Centroid(1)  , x_d(:,2)-Centroid(2)]; Dr2=Dr(:,1).^2+Dr(:,2).^2;
   I_t=sum(Dr2)/length(x_d(:,1));
   I(it)=I_t/I0;
   
   
   figure(4);
   plot(Dt*[1:it], I(1:it), 'r'); 
   title('Particle dispersion');
   xlabel('t'); ylabel('I/I_0'); drawnow;

   %%%%%%%%%%%%%%%%%%%%%%%%%
% Output
    if mod(it,100) == 0
        t = it*Dt;
% Verifica sulla divergenza alla fine dello step
        MaxDiv = max(max(abs(DivCalc(U,V))));
       
        disp(['Massima divergenza sul campo = ',num2str(MaxDiv)])
% Interpolazione nei corners della velocità
        i = 1:Nx;     j = 2:Ny-1;   Iu(i,j)  = (U(i,j) + U(i,j-1))/2; 
        Iu(i,1) = Udwn;             Iu(i,Ny) = Uup;
        i = 2:Nx-1;   j = 1:Ny;     Iv(i,j)  = (V(i,j) + V(i-1,j))/2; 
        Iv(1,j) = VWest;            Iv(Nx,j) = VEst;
% Mesh 2D
        X = repmat(x,Ny,1);         Y = repmat(y',1,Nx);
% Grafica
        figure(2); clf; 
        streamslice(X,Y,Iu',Iv');     axis([0 Lx 0 Ly ]);     axis square;   
        title(['Harlow-Welch. Driven cavity. t = ',num2str(t)]); drawnow;
        figure(3); clf;
        PSI  = GivePsi(U,V);                 % Calcolo della Psi da U e V
        vneg = linspace(min(min(PSI)),0,20);
        vpos = linspace(0,max(max(PSI)),20);
        contour(x,y,PSI',vneg,'k'); axis square; hold on;
        contour(x,y,PSI',vpos,'r'); 
        title('Harlow-Welch. Driven cavity. Streamlines da \psi');
        drawnow; hold off
    end
    figure (5)
    if it == round(Nt/3) 
        
        subplot(1,3,1)
        vneg = linspace(min(min(PSI)),0,20);
        vpos = linspace(0,max(max(PSI)),20);
        contour(x,y,PSI',vneg,'k'); axis square; hold on;
        contour(x,y,PSI',vpos,'r'); 
        title(['Harlow-Welch. Driven cavity. t = ',num2str(t)])
        drawnow; hold off

    elseif it == round(2*Nt/3)
        subplot(1,3,2)
        vneg = linspace(min(min(PSI)),0,20);
        vpos = linspace(0,max(max(PSI)),20);
        contour(x,y,PSI',vneg,'k'); axis square; hold on;
        contour(x,y,PSI',vpos,'r'); 
        title(['Harlow-Welch. Driven cavity. t = ',num2str(t)])
        drawnow; hold off

     elseif it == round(Nt)
        subplot(1,3,3)
        vneg = linspace(min(min(PSI)),0,20);
        vpos = linspace(0,max(max(PSI)),20);
        contour(x,y,PSI',vneg,'k'); axis square; hold on;
        contour(x,y,PSI',vpos,'r'); 
        title(['Harlow-Welch. Driven cavity. t = ',num2str(t)])
        drawnow; hold off

        
    end


        figure (6)
    if it == round(Nt/3) 
        subplot(1,3,1)
   s1=scatter(x_d(:,1), x_d(:,2), 20, ColorDots ,'filled');
   hold on
   axis([0 Lx 0 Ly]); axis square;
   title(['Particle scattering t= ', num2str(it*Dt)]); drawnow;

    elseif it == round(2*Nt/3)
        subplot(1,3,2)
   s1=scatter(x_d(:,1), x_d(:,2), 20, ColorDots ,'filled');
   hold on
   axis([0 Lx 0 Ly]); axis square;
   title(['Particle scattering t= ', num2str(it*Dt)]); drawnow;

     elseif it == round(Nt)
         subplot(1,3,3)
   s1=scatter(x_d(:,1), x_d(:,2), 20, ColorDots ,'filled');
   hold on
   axis([0 Lx 0 Ly]); axis square;
   title(['Particle scattering t= ', num2str(it*Dt)]); drawnow;

        
    end
    
end


