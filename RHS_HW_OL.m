function [Fu,Fv] = RHS_HW_OL(U,V,Uup, Udwn)

global Nx Ny Re h hq UWest UEst VNord VSud VWest VEst

Fu = zeros(size(U)); Fv = zeros(size(V));
Du = zeros(size(U)); Dv = zeros(size(V));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convective term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U COMPONENT

% dx(Ix(u)*Ix(u))
i = 1:Nx-1;  j = 1:Ny-1;
Iu(i,j) = (U(i+1,j) + U(i,j))/2;              %cell centered mean value
Iuq     = Iu.*Iu;                             
i = 2:Nx-1; j = 1:Ny-1;
Fu(i,j) = Fu(i,j) + (Iuq(i,j)-Iuq(i-1,j))/h;  

% dy(Iy(u)*Ix(v))
Iu = zeros(Nx,Ny); Iv = zeros(Nx,Ny);
Iu(:,Ny) = Uup;  Iu(:,1)  = Udwn; Iu(1,:)  = UWest;  Iu(Nx,:) = UEst; %% Uup %% Udwn
Iv(:,Ny) = VNord;  Iv(:,1)  = VSud; Iv(1,:)  = VWest;  Iv(Nx,:) = VEst;
i = 2:Nx-1; j = 2:Ny-1;
Iu(i,j) = (U(i,j) + U(i,j-1))/2;              % mean value on the corners
Iv(i,j) = (V(i,j) + V(i-1,j))/2;
IuIv = Iu.*Iv;                                
i = 2:Nx-1; j = 1:Ny-1;
Fu(i,j) = Fu(i,j) + (IuIv(i,j+1)-IuIv(i,j))/h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V COMPONENT

% dx(Iy(u)*Ix(v))
i = 1:Nx-1; j = 2:Ny-1;
Fv(i,j) = Fv(i,j) + (IuIv(i+1,j)-IuIv(i,j))/h;

% dy(Iy(v)*Iy(v))
i = 1:Nx-1; j = 1:Ny-1;
Iv(i,j) = (V(i,j+1) + V(i,j))/2;              
Ivq = Iv.*Iv;                                 
i = 1:Nx-1; j = 2:Ny-1;
Fv(i,j) = Fv(i,j) + (Ivq(i,j)-Ivq(i,j-1))/h;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusive Term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U COMPONENT

i = 2:Nx-1; j = 1;
Du(i,j) = (U(i,j+1) + U(i+1,j) + U(i-1,j) + 2*Udwn  - 5*U(i,j))/hq; %% Udwn
i = 2:Nx-1; j = Ny-1;
Du(i,j) = (U(i,j-1) + U(i+1,j) + U(i-1,j) + 2*Uup - 5*U(i,j))/hq;   %% Uup
i = 2:Nx-1; j = 2:Ny-2;
Du(i,j) = (U(i,j-1) + U(i,j+1) + U(i+1,j) + U(i-1,j) - 4*U(i,j))/hq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V COMPONENT
i = 1;    j = 2:Ny-1;
Dv(i,j) = (V(i,j+1) + V(i+1,j) + V(i,j-1) + 2*VWest  - 5*V(i,j))/hq;
i = Nx-1; j = 2:Ny-1;
Dv(i,j) = (V(i,j+1) + V(i-1,j) + V(i,j-1) + 2*VEst   - 5*V(i,j))/hq;
i = 2:Nx-2; j = 2:Ny-1;
Dv(i,j) = (V(i,j-1) + V(i,j+1) + V(i+1,j) + V(i-1,j) - 4*V(i,j))/hq;


Fu = -Fu + Du/Re;
Fv = -Fv + Dv/Re;

end

