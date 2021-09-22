clear all; clc;

global S0 K r T vol Smax Ntau dS_1 dS_2 dS_3 corr face q dummy barrier grid ki_grid 

vol = 0.3; vol_x=vol;vol_y=vol;vol_z=vol;
corr = 0.5; corr_xy=corr;corr_yz=corr;corr_zx=corr;

S0=100; K=[85, 85, 90, 90, 95, 95]; r = 0.03; T = 3;  Smax=150; 
dS_1= 15; dS_2= 5; dS_3= 20; corr = 0.5; face = 100;
q = 0.3:-0.05:0.05; dummy = 0.26; barrier = 50; 

Sx = get_non_uniform_3d(Smax,dS_1,dS_2,dS_3); Sy = get_non_uniform_3d(Smax,dS_1,dS_2,dS_3); Sz = get_non_uniform_3d(Smax,dS_1,dS_2,dS_3);
dSx= diff(Sx); dSy= diff(Sy); dSz= diff(Sz);
Nx= numel(Sx); Ny= numel(Sy); Nz = numel(Sz);
dtau =  1/30; Ntau = T/dtau;
step = Ntau * (T + 1);  
grid = zeros(Nx,Ny,Nz); %knock-in barrier로 떨어지지 않는 경우 옵션가격
ki_grid = zeros(Nx,Ny,Nz);




%initial_conditions = pay_off
for idx_x= 1:Nx   
for idx_y= 1:Ny
for idx_z= 1:Nz
if (Sx(idx_x) <= barrier|| Sx(idx_y) <= barrier || Sx(idx_z) <= barrier)
    grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
    ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
                         
elseif (Sx(idx_x) < K(1)||Sy(idx_y) < K(1)||Sz(idx_z) < K(1))
       grid(idx_x, idx_y, idx_z) = face * (1 + dummy);
       ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
else
       grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       ki_grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       
end
end
end
end

%traverse time
tag=1;
for iter=1:Ntau
    
if iter == step(tag)
    gx=min(find(x >= x0*K(tag+1)));
    gy=min(find(y >= y0*K(tag+1)));
    gz=min(find(z >= z0*K(tag+1)));
    grid(gx: Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    ki_grid(gx:Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    
    tag= tag+1;
end


%update_ki_grid
for i=2:Nx - 1
for j=2:Ny - 1
for k=2:Nz - 1
S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * ki_grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * ki_grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * ki_grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * ki_grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(ki_grid(i+1,j+1,k) - ki_grid(i-1,j+1,k) - ki_grid(i+1,j-1,k) + ki_grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(ki_grid(i,j+1,k+1) - ki_grid(i,j-1,k+1) - ki_grid(i,j+1,k-1) + ki_grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(ki_grid(i+1,j,k+1) - ki_grid(i+1,j,k-1)- ki_grid(i-1,j,k+1) + ki_grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * ki_grid(i,j,k)


ki_grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_ki_grid
for i=3:Nx-1
for j=3:Ny-1
for k=3:Nz-1
    ki_grid(Nx-1, 1:Ny-1, 1:Nz-1) = ki_grid(Nx-2, 1:Ny-1, 1:Nz-1);
    ki_grid(1:Nx-1, Ny-1, 1:Nz-1) = ki_grid(1:Nx-1, Ny-2, 1:Nz-1);
    ki_grid(1:Nx-1, 1:Ny-1, Nz-1) = ki_grid(1:Nx-1 , 1:Ny-1, Nz-2);
end;end;end


%traverse grid
minimum_idx_x = get_argmin(Sx, barrier);
minimum_idx_y = get_argmin(Sy, barrier);
minimum_idx_z = get_argmin(Sz, barrier);
grid(minimum_idx_x, :, :) = ki_grid(minimum_idx_x, :, :);
grid(:, minimum_idx_y, :) = ki_grid(:, minimum_idx_y, :);
grid(:, :, minimum_idx_z) = ki_grid(:, :, minimum_idx_z);

%update_grid

for i=2:Nx - 1
for j=2:Ny - 1
for k=2:Nz - 1
S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(grid(i+1,j+1,k) - grid(i-1,j+1,k) - grid(i+1,j-1,k) + grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(grid(i,j+1,k+1) - grid(i,j-1,k+1) - grid(i,j+1,k-1) + grid(i,j-1,k-1));
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(grid(i+1,j,k+1) - grid(i+1,j,k-1)- grid(i-1,j,k+1) + grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * grid(i,j,k)


grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_grid
for i=3:Nx-1
for j=3:Ny-1
for k=3:Nz-1
    grid(Nx-1, 1:Ny-1, 1:Nz-1) = grid(Nx-2, 1:Ny-1, 1:Nz-1);
    grid(1:Nx-1, Ny-1, 1:Nz-1) = grid(1:Nx-1, Ny-2, 1:Nz-1);
    grid(1:Nx-1, 1:Ny-1, Nz-1) = grid(1:Nx-1 , 1:Ny-1, Nz-2);
end;end;end
end



%coefficients_i_k_j

%initial_conditions = pay_off
for idx_x= 1:Nx   
for idx_z= 1:Nz
for idx_y= 1:Ny
if (Sx(idx_x) <= barrier|| Sx(idx_y) <= barrier || Sx(idx_z) <= barrier)
    grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
    ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
                         
elseif (Sx(idx_x) < K(1)||Sy(idx_y) < K(1)||Sz(idx_z) < K(1))
       grid(idx_x, idx_y, idx_z) = face * (1 + dummy);
       ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
else
       grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       ki_grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       
end
end
end
end

%traverse time
tag=1;
for iter=1:Ntau
    
if iter == step(tag)
    gx=min(find(x >= x0*K(tag+1)));
    gy=min(find(y >= y0*K(tag+1)));
    gz=min(find(z >= z0*K(tag+1)));
    grid(gx: Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    ki_grid(gx:Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    
    tag= tag+1;
end


%update_ki_grid

for i=2:Nx - 1
for k=2:Nz - 1
for j=2:Ny - 1

S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * ki_grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * ki_grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * ki_grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * ki_grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(ki_grid(i+1,j+1,k) - ki_grid(i-1,j+1,k) - ki_grid(i+1,j-1,k) + ki_grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(ki_grid(i,j+1,k+1) - ki_grid(i,j-1,k+1) - ki_grid(i,j+1,k-1) + ki_grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(ki_grid(i+1,j,k+1) - ki_grid(i+1,j,k-1)- ki_grid(i-1,j,k+1) + ki_grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * ki_grid(i,j,k)


ki_grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_ki_grid
for i=3:Nx-1
for k=3:Nz-1
for j=3:Ny-1
    ki_grid(Nx-1,1:Ny-1, 1:Nz-1) = ki_grid(Nx-2, 1:Ny-1, 1:Nz-1);
    ki_grid(1:Nx-1, Ny-1, 1:Nz-1) = ki_grid(1:Nx-1,Ny-2, 1:Nz-1);
    ki_grid(1:Nx-1, 1:Ny-1, Nz-1) = ki_grid(1:Nx-1,1:Ny-1, Nz-2);
end;end;end


%traverse grid
minimum_idx_x = get_argmin(Sx, barrier);
minimum_idx_y = get_argmin(Sy, barrier);
minimum_idx_z = get_argmin(Sz, barrier);
grid(minimum_idx_x, :, :) = ki_grid(minimum_idx_x, :, :);
grid(:, minimum_idx_y, :) = ki_grid(:, minimum_idx_y, :);
grid(:, :, minimum_idx_z) = ki_grid(:, :, minimum_idx_z);

%update_grid

for i=2:Nx - 1
for k=2:Nz - 1
for j=2:Ny - 1

S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(grid(i+1,j+1,k) - grid(i-1,j+1,k) - grid(i+1,j-1,k) + grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(grid(i,j+1,k+1) - grid(i,j-1,k+1) - grid(i,j+1,k-1) + grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(grid(i+1,j,k+1) - grid(i+1,j,k-1)- grid(i-1,j,k+1) + grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * grid(i,j,k)


grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_grid
for i=3:Nx-1
for k=3:Nz-1
for j=3:Ny-1
    grid(Nx-1, 1:Ny-1, 1:Nz-1) = grid(Nx-2, 1:Ny-1, 1:Nz-1);
    grid(1:Nx-1, Ny-1, 1:Nz-1) = grid(1:Nx-1, Ny-2, 1:Nz-1);
    grid(1:Nx-1, 1:Ny-1, Nz-1) = grid(1:Nx-1 , 1:Ny-1, Nz-2);
end;end;end

end


%coefficients_k_i_j

%initial_conditions = pay_off
for idx_z= 1:Nz
for idx_x= 1:Nx   
for idx_y= 1:Ny
if (Sx(idx_x) <= barrier|| Sx(idx_y) <= barrier || Sx(idx_z) <= barrier)
    grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
    ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
                         
elseif (Sx(idx_x) < K(1)||Sy(idx_y) < K(1)||Sz(idx_z) < K(1))
       grid(idx_x, idx_y, idx_z) = face * (1 + dummy);
       ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
else
       grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       ki_grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       
end
end
end
end

%traverse time
tag=1;
for iter=1:Ntau
    
if iter == step(tag)
    gx=min(find(x >= x0*K(tag+1)));
    gy=min(find(y >= y0*K(tag+1)));
    gz=min(find(z >= z0*K(tag+1)));
    grid(gx: Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    ki_grid(gx:Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    
    tag= tag+1;
end


%update_ki_grid

for k=2:Nz - 1
for i=2:Nx - 1
for j=2:Ny - 1

S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * ki_grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * ki_grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * ki_grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * ki_grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(ki_grid(i+1,j+1,k) - ki_grid(i-1,j+1,k) - ki_grid(i+1,j-1,k) + ki_grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(ki_grid(i,j+1,k+1) - ki_grid(i,j-1,k+1) - ki_grid(i,j+1,k-1) + ki_grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(ki_grid(i+1,j,k+1) - ki_grid(i+1,j,k-1)- ki_grid(i-1,j,k+1) + ki_grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * ki_grid(i,j,k)


ki_grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_ki_grid
for k=3:Nz-1
for i=3:Nx-1
for j=3:Ny-1
    ki_grid(Nx-1,1:Ny-1, 1:Nz-1) = ki_grid(Nx-2, 1:Ny-1, 1:Nz-1);
    ki_grid(1:Nx-1, Ny-1, 1:Nz-1) = ki_grid(1:Nx-1,Ny-2, 1:Nz-1);
    ki_grid(1:Nx-1, 1:Ny-1, Nz-1) = ki_grid(1:Nx-1,1:Ny-1, Nz-2);
end;end;end


%traverse grid
minimum_idx_x = get_argmin(Sx, barrier);
minimum_idx_y = get_argmin(Sy, barrier);
minimum_idx_z = get_argmin(Sz, barrier);
grid(minimum_idx_x, :, :) = ki_grid(minimum_idx_x, :, :);
grid(:, minimum_idx_y, :) = ki_grid(:, minimum_idx_y, :);
grid(:, :, minimum_idx_z) = ki_grid(:, :, minimum_idx_z);

%update_grid
for k=2:Nz - 1
for i=2:Nx - 1
for j=2:Ny - 1

S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(grid(i+1,j+1,k) - grid(i-1,j+1,k) - grid(i+1,j-1,k) + grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(grid(i,j+1,k+1) - grid(i,j-1,k+1) - grid(i,j+1,k-1) + grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(grid(i+1,j,k+1) - grid(i+1,j,k-1)- grid(i-1,j,k+1) + grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * grid(i,j,k)


grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_grid
for k=3:Nz-1
for i=3:Nx-1
for j=3:Ny-1

    grid(Nx-1, 1:Ny-1, 1:Nz-1) = grid(Nx-2, 1:Ny-1, 1:Nz-1);
    grid(1:Nx-1, Ny-1, 1:Nz-1) = grid(1:Nx-1, Ny-2, 1:Nz-1);
    grid(1:Nx-1, 1:Ny-1, Nz-1) = grid(1:Nx-1 , 1:Ny-1, Nz-2);
end;end;end

end


%coefficients_i_k_j
%initial_conditions = pay_off
for idx_x= 1:Nx   
for idx_z= 1:Nz
for idx_y= 1:Ny
if (Sx(idx_x) <= barrier|| Sx(idx_y) <= barrier || Sx(idx_z) <= barrier)
    grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
    ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
                         
elseif (Sx(idx_x) < K(1)||Sy(idx_y) < K(1)||Sz(idx_z) < K(1))
       grid(idx_x, idx_y, idx_z) = face * (1 + dummy);
       ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
else
       grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       ki_grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       
end
end
end
end

%traverse time
tag=1;
for iter=1:Ntau
    
if iter == step(tag)
    gx=min(find(x >= x0*K(tag+1)));
    gy=min(find(y >= y0*K(tag+1)));
    gz=min(find(z >= z0*K(tag+1)));
    grid(gx: Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    ki_grid(gx:Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    
    tag= tag+1;
end


%update_ki_grid

for i=2:Nx - 1
for k=2:Nz - 1
for j=2:Ny - 1

S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * ki_grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * ki_grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * ki_grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * ki_grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(ki_grid(i+1,j+1,k) - ki_grid(i-1,j+1,k) - ki_grid(i+1,j-1,k) + ki_grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(ki_grid(i,j+1,k+1) - ki_grid(i,j-1,k+1) - ki_grid(i,j+1,k-1) + ki_grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(ki_grid(i+1,j,k+1) - ki_grid(i+1,j,k-1)- ki_grid(i-1,j,k+1) + ki_grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * ki_grid(i,j,k)


ki_grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_ki_grid
for i=3:Nx-1
for k=3:Nz-1
for j=3:Ny-1
    ki_grid(Nx-1,1:Ny-1, 1:Nz-1) = ki_grid(Nx-2, 1:Ny-1, 1:Nz-1);
    ki_grid(1:Nx-1, Ny-1, 1:Nz-1) = ki_grid(1:Nx-1,Ny-2, 1:Nz-1);
    ki_grid(1:Nx-1, 1:Ny-1, Nz-1) = ki_grid(1:Nx-1,1:Ny-1, Nz-2);
end;end;end


%traverse grid
minimum_idx_x = get_argmin(Sx, barrier);
minimum_idx_y = get_argmin(Sy, barrier);
minimum_idx_z = get_argmin(Sz, barrier);
grid(minimum_idx_x, :, :) = ki_grid(minimum_idx_x, :, :);
grid(:, minimum_idx_y, :) = ki_grid(:, minimum_idx_y, :);
grid(:, :, minimum_idx_z) = ki_grid(:, :, minimum_idx_z);

%update_grid

for i=2:Nx - 1
for k=2:Nz - 1
for j=2:Ny - 1

S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(grid(i+1,j+1,k) - grid(i-1,j+1,k) - grid(i+1,j-1,k) + grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(grid(i,j+1,k+1) - grid(i,j-1,k+1) - grid(i,j+1,k-1) + grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(grid(i+1,j,k+1) - grid(i+1,j,k-1)- grid(i-1,j,k+1) + grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * grid(i,j,k)


grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_grid
for i=3:Nx-1
for k=3:Nz-1
for j=3:Ny-1

    grid(Nx-1, 1:Ny-1, 1:Nz-1) = grid(Nx-2, 1:Ny-1, 1:Nz-1);
    grid(1:Nx-1, Ny-1, 1:Nz-1) = grid(1:Nx-1, Ny-2, 1:Nz-1);
    grid(1:Nx-1, 1:Ny-1, Nz-1) = grid(1:Nx-1 , 1:Ny-1, Nz-2);
end;end;end

end

%coefficients_j_i_k

%initial_conditions = pay_off
for idx_y= 1:Ny
for idx_x= 1:Nx   
for idx_z= 1:Nz

if (Sx(idx_x) <= barrier|| Sx(idx_y) <= barrier || Sx(idx_z) <= barrier)
    grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
    ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
                         
elseif (Sx(idx_x) < K(1)||Sy(idx_y) < K(1)||Sz(idx_z) < K(1))
       grid(idx_x, idx_y, idx_z) = face * (1 + dummy);
       ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
else
       grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       ki_grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       
end
end
end
end

%traverse time
tag=1;
for iter=1:Ntau
    
if iter == step(tag)
    gx=min(find(x >= x0*K(tag+1)));
    gy=min(find(y >= y0*K(tag+1)));
    gz=min(find(z >= z0*K(tag+1)));
    grid(gx: Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    ki_grid(gx:Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    
    tag= tag+1;
end


%update_ki_grid

for i=2:Nx - 1
for j=2:Ny - 1
for k=2:Nz - 1
S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * ki_grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * ki_grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * ki_grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * ki_grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(ki_grid(i+1,j+1,k) - ki_grid(i-1,j+1,k) - ki_grid(i+1,j-1,k) + ki_grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(ki_grid(i,j+1,k+1) - ki_grid(i,j-1,k+1) - ki_grid(i,j+1,k-1) + ki_grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(ki_grid(i+1,j,k+1) - ki_grid(i+1,j,k-1)- ki_grid(i-1,j,k+1) + ki_grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * ki_grid(i,j,k)


ki_grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_ki_grid
for j=3:Ny-1
for i=3:Nx-1
for k=3:Nz-1

    ki_grid(Nx-1,1:Ny-1, 1:Nz-1) = ki_grid(Nx-2, 1:Ny-1, 1:Nz-1);
    ki_grid(1:Nx-1, Ny-1, 1:Nz-1) = ki_grid(1:Nx-1,Ny-2, 1:Nz-1);
    ki_grid(1:Nx-1, 1:Ny-1, Nz-1) = ki_grid(1:Nx-1,1:Ny-1, Nz-2);
end;end;end

%traverse grid
minimum_idx_x = get_argmin(Sx, barrier);
minimum_idx_y = get_argmin(Sy, barrier);
minimum_idx_z = get_argmin(Sz, barrier);
grid(minimum_idx_x, :, :) = ki_grid(minimum_idx_x, :, :);
grid(:, minimum_idx_y, :) = ki_grid(:, minimum_idx_y, :);
grid(:, :, minimum_idx_z) = ki_grid(:, :, minimum_idx_z);

%update_grid
for j=2:Ny - 1
for i=2:Nx - 1
for k=2:Nz - 1
S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(grid(i+1,j+1,k) - grid(i-1,j+1,k) - grid(i+1,j-1,k) + grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(grid(i,j+1,k+1) - grid(i,j-1,k+1) - grid(i,j+1,k-1) + grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(grid(i+1,j,k+1) - grid(i+1,j,k-1)- grid(i-1,j,k+1) + grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * grid(i,j,k)


grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_grid
for j=3:Ny-1
for i=3:Nx-1
for k=3:Nz-1
    grid(Nx-1, 1:Ny-1, 1:Nz-1) = grid(Nx-2, 1:Ny-1, 1:Nz-1);
    grid(1:Nx-1, Ny-1, 1:Nz-1) = grid(1:Nx-1, Ny-2, 1:Nz-1);
    grid(1:Nx-1, 1:Ny-1, Nz-1) = grid(1:Nx-1 , 1:Ny-1, Nz-2);
end;end;end
end


%coefficients_k_j_i
%initial_conditions = pay_off
   
for idx_z= 1:Nz
for idx_y= 1:Ny
for idx_x= 1:Nx
if (Sx(idx_x) <= barrier|| Sx(idx_y) <= barrier || Sx(idx_z) <= barrier)
    grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
    ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
                         
elseif (Sx(idx_x) < K(1)||Sy(idx_y) < K(1)||Sz(idx_z) < K(1))
       grid(idx_x, idx_y, idx_z) = face * (1 + dummy);
       ki_grid(idx_x, idx_y, idx_z) = face * min([Sx(idx_x), Sy(idx_y), Sz(idx_z)]) / S0;
else
       grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       ki_grid(idx_x, idx_y, idx_z) = face * (1 + q(1));
       
end
end
end
end

%traverse time
tag=1;
for iter=1:Ntau
    
if iter == step(tag)
    gx=min(find(x >= x0*K(tag+1)));
    gy=min(find(y >= y0*K(tag+1)));
    gz=min(find(z >= z0*K(tag+1)));
    grid(gx: Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    ki_grid(gx:Nx, gy:Ny, gz:Nz) = face*(1+q(tag+1));
    
    tag= tag+1;
end


%update_ki_grid
for k=2:Nz - 1
for j=2:Ny - 1    
for i=2:Nx - 1
   
S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * ki_grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * ki_grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * ki_grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * ki_grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * ki_grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(ki_grid(i+1,j+1,k) - ki_grid(i-1,j+1,k) - ki_grid(i+1,j-1,k) + ki_grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(ki_grid(i,j+1,k+1) - ki_grid(i,j-1,k+1) - ki_grid(i,j+1,k-1) + ki_grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(ki_grid(i+1,j,k+1) - ki_grid(i+1,j,k-1)- ki_grid(i-1,j,k+1) + ki_grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * ki_grid(i,j,k)


ki_grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_ki_grid

for k=3:Nz-1
for j=3:Ny-1
for i=3:Nx-1
    ki_grid(Nx-1,1:Ny-1, 1:Nz-1) = ki_grid(Nx-2, 1:Ny-1, 1:Nz-1);
    ki_grid(1:Nx-1, Ny-1, 1:Nz-1) = ki_grid(1:Nx-1,Ny-2, 1:Nz-1);
    ki_grid(1:Nx-1, 1:Ny-1, Nz-1) = ki_grid(1:Nx-1,1:Ny-1, Nz-2);
end;end;end


%traverse grid
minimum_idx_x = get_argmin(Sx, barrier);
minimum_idx_y = get_argmin(Sy, barrier);
minimum_idx_z = get_argmin(Sz, barrier);
grid(minimum_idx_x, :, :) = ki_grid(minimum_idx_x, :, :);
grid(:, minimum_idx_y, :) = ki_grid(:, minimum_idx_y, :);
grid(:, :, minimum_idx_z) = ki_grid(:, :, minimum_idx_z);

%update_grid
for k=2:Nz - 1
for j=2:Ny - 1
for i=2:Nx - 1    
    
S(1)= 1 / dtau + 0.5 * r + (vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i)*(dSx(i-1) + dSx(i))) +  (vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j)*(dSy(j-1) + dSy(j)))...
    +  (vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k)*(dSz(k-1) + dSz(k)));
S(2)= ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) * grid(i+1,j,k);
S(3)= ((vol_x^2 * Sx(i)^2 - r * Sx(i) * dSx(i)) / (dSx(i-1)*(dSx(i-1) + dSx(i)))) * grid(i,j+1,k);
S(4)= ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) * grid(i,j,k+1);
S(5)= ((vol_y^2 * Sy(j)^2 - r * Sy(j) * dSy(j)) / (dSy(j-1)*(dSy(j-1) + dSy(j)))) * grid(i-1,j,k);
S(6)= ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k)))) * grid(i,j-1,k);
S(7)= ((vol_z^2 * Sz(k)^2 - r * Sz(k) * dSz(k)) / (dSz(k-1)*(dSz(k-1) + dSz(k)))) * grid(i,j,k-1);
S(8)= (corr_xy * vol_x * vol_y * Sx(i) * Sy(j) / (dSx(i) * dSy(j) + dSx(i-1) * dSy(j) + dSx(i) * dSy(j-1) + dSx(i-1) * dSy(j-1)))...
    *(grid(i+1,j+1,k) - grid(i-1,j+1,k) - grid(i+1,j-1,k) + grid(i-1,j-1,k))
S(9)= (corr_yz * vol_y * vol_z * Sy(j) * Sz(k) / (dSy(j) * dSz(k) + dSy(j-1) * dSz(k) + dSy(j) * dSz(k-1) + dSy(j-1) * dSz(k-1)))...
    *(grid(i,j+1,k+1) - grid(i,j-1,k+1) - grid(i,j+1,k-1) + grid(i,j-1,k-1))
S(10)= (corr_zx * vol_z * vol_x * Sz(k) * Sx(i) / (dSz(k) * dSx(i) + dSz(k-1) * dSx(i) + dSz(k) * dSx(i-1) + dSz(k-1) * dSx(i-1)))...
    *(grid(i+1,j,k+1) - grid(i+1,j,k-1)- grid(i-1,j,k+1) + grid(i-1,j,k-1))
S(11)= (((1 / dtau - 0.5*r)  - ((vol_x^2 * Sx(i)^2 + r * Sx(i) * dSx(i-1)) / (dSx(i)*(dSx(i-1) + dSx(i)))) - ((vol_y^2 * Sy(j)^2 + r * Sy(j) * dSy(j-1)) / (dSy(j)*(dSy(j-1) + dSy(j)))) - ((vol_z^2 * Sz(k)^2 + r * Sz(k) * dSz(k-1)) / (dSz(k)*(dSz(k-1) + dSz(k))))...
    - r * Sx(i) * dSx(i-1) - r * Sy(j) * dSy(j-1) - r * Sz(k) * dSz(k-1))) * grid(i,j,k)


grid_ijk= (S(2)+S(3)+S(4)+S(5)+S(6)+S(7)+S(8)+S(9)+S(10)+S(11))/S(1)
end;end;end


%neumann_boundary_grid 
for k=3:Nz-1
for j=3:Ny-1
for i=3:Nx-1
    grid(Nx-1, 1:Ny-1, 1:Nz-1) = grid(Nx-2, 1:Ny-1, 1:Nz-1);
    grid(1:Nx-1, Ny-1, 1:Nz-1) = grid(1:Nx-1, Ny-2, 1:Nz-1);
    grid(1:Nx-1, 1:Ny-1, Nz-1) = grid(1:Nx-1 , 1:Ny-1, Nz-2);
end;end;end
end


%price
% price = interp3(i,j,k,grid,Nx,Ny,Nz)

%non_uniform_condition
function domain = get_non_uniform_3d(Smax,dS_1,dS_2,dS_3)
domain_1 = 0:dS_1:round(Smax/3);
domain_2 = round(Smax/3):dS_2:round(Smax*2/3);
domain_3 = round(Smax*2/3):dS_3:round(Smax);
domain= unique(cat(2,domain_1,domain_2,domain_3));
end 


% argmin 
function minimum_idx = get_argmin(arr, limit)
valid_idx = find(arr >= limit);
minimum_idx = find(valid_idx==min(valid_idx));
end

    





    



 
