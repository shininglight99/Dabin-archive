clear;clc;

K= 230; Smax= 800; vol_x= 0.5; vol_y=vol_x; r=0.03; T=0.01; dt=1/10000; Nx=100; Ny=Nx; Nt=T/dt; rho=0.5;

Sx = linspace(0,Smax,Nx);Sy=Sx; h_x = diff(Sx);h_y=h_x;


u(1:Nx,1:Ny,1:Nt+1)=0; 
for i=1:Nx
    for j=1:Ny
        u(i,j,1) = max(min(Sx(i),Sy(j))-K,0); 
    end
end
figure(1);clf;hold on; mesh(Sx,Sy,u(:,:,1)');view(3);



%set coefficient
for i=2:Nx-1
    for j=2:Ny-1
        A = ((vol_x*Sx(i)).^2 - r*Sx(i).*h_x(i))./(h_x(i-1).*(h_x(i-1)+h_x(i)))...
            + ((vol_y*Sy(j)).^2 - r*Sy(j).*h_y(j))./(h_y(j-1).*(h_y(j-1)+h_y(j)))+ 0.5*r + 1/dt;
        B = ((vol_x*Sx(i)).^2 + r*Sx(i).*h_x(i-1))./(h_x(i).*(h_x(i-1)+h_x(i)));
        C = ((vol_y*Sy(j)).^2 + r*Sy(j).*h_y(j-1))./(h_y(j).*(h_y(j-1)+h_y(j)));
        D = ((vol_x*Sx(i)).^2)./(h_x(i-1).*(h_x(i-1)+h_x(i)));
        E = ((vol_y*Sy(j)).^2)./(h_y(j-1).*(h_y(j-1)+h_y(j)));
        F = (- r*Sx(i).*h_x(i))./(h_x(i-1).*(h_x(i-1)+h_x(i)));
        G = (- r*Sy(j).*h_y(j))./(h_y(j-1).*(h_y(j-1)+h_y(j)));
        H = -((vol_x*Sx(i)).^2 + r*Sx(i).*h_x(i-1))./(h_x(i).*(h_x(i-1)+h_x(i)))...
            -((vol_y*Sy(j)).^2 + r*Sy(j).*h_y(j-1))./(h_y(j).*(h_y(j-1)+h_y(j)))- 0.5*r + 1/dt;
        I = (rho*vol_x*vol_y*Sx(i).*Sy(j))./((h_x(i).*h_y(j)) + (h_x(i-1).*h_y(j))...
            +(h_x(i).*h_y(j-1))+(h_x(i-1).*h_y(j-1)));
    end
end

B_ = B./A; C_ = C./A; D_ = D./A; E_ = E./A; F_ = F./A; 
G_ = G./A; H_ = H./A; I_ = I./A;

%loop

for n=1:Nt
    for i=2:Nx-1
        for j=2:Ny-1
            u(i,j,n+1) = B_*u(i+1,j,n) + C_*u(i,j+1,n) + D_*u(i-1,j,n) + E_*u(i,j-1,n)...
                - F_*u(i-1,j,n+1) - G_*u(i,j-1,n+1) + H_*u(i,j,n) + I_*(u(i+1,j+1,n)-u(i-1,j+1,n)-u(i+1,j-1,n)+u(i-1,j-1,n));
        
        end
    end
    u(2:Nx-1,Ny,n+1)=2*u(2:Nx-1,Ny-1,n+1)-u(2:Nx-1,Ny-2,n+1);
    u(Nx,2:Ny,n+1)=2*u(Nx-1,2:Ny,n+1)-u(Nx-2,2:Ny,n+1);
    
    for j=2:Ny-1
        for i=2:Nx-1
            u(i,j,n+1) = B_*u(i+1,j,n) + C_*u(i,j+1,n) + D_*u(i-1,j,n) + E_*u(i,j-1,n)...
                - F_*u(i-1,j,n+1) - G_*u(i,j-1,n+1) + H_*u(i,j,n) + I_*(u(i+1,j+1,n)-u(i-1,j+1,n)-u(i+1,j-1,n)+u(i-1,j-1,n));
        
        end
    end
    u(2:Nx-1,Ny,n+1)=2*u(2:Nx-1,Ny-1,n+1)-u(2:Nx-1,Ny-2,n+1);
    u(Nx,2:Ny,n+1)=2*u(Nx-1,2:Ny,n+1)-u(Nx-2,2:Ny,n+1);
    
   
    figure(1);clf; mesh(Sx,Sy,u(:,:,n+1)');
    
end

