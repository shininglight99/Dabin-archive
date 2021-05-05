#!/usr/bin/env python
# coding: utf-8

# In[16]:


import numpy as np
from itertools import permutations
import copy
import time
import timeit
from utils import _get_non_uniform_3D_S, _get_argmin,  _printProgress #,_update_grid


class FDM3D:
    
    def __init__(self, S0, K, r, T, sigma, Smax, Ntau, dS, corr, face, q, dummy, barrier):
        self.S0x, self.S0y, self.S0z = S0
        self.Kx, self.Ky, self.Kz = K  # K = ([85, 85, 90, 90, 95, 95], [85, 85, 90, 90, 95, 95],[85, 85, 90, 90, 95, 95])
        self.r = r   # r= 0.03
        self.T = T   # T = 3
        self.sigma_x, self.sigma_y, self.sigma_z = sigma   # [0.3, 0.3, 0.3]
        self.Smax_x, self.Smax_y, self.Smax_z = Smax   # [150, 150, 150]
        self.dS = dS   # dS in three sub-domain
        self.Ntau = Ntau   # T x 360
        self.corr_xy, self.corr_yz, self.corr_zx = corr
        self.face = face
        self.q = q   # coupon rate
        self.dummy = dummy   # dummy rate
        self.barrier = barrier   # knock-in barrier
        
        self.Sx, self.Sy, self.Sz  = _get_non_uniform_3D_S(self.Smax_x, self.dS), _get_non_uniform_3D_S(self.Smax_y, self.dS), _get_non_uniform_3D_S(self.Smax_z, self.dS)
        
        
        self.dSx, self.dSy, self.dSz = self.Sx[1:] - self.Sx[:-1], self.Sy[1:] - self.Sy[:-1], self.Sz[1:] - self.Sz[:-1]
        self.Nx, self.Ny, self.Nz = len(self.Sx), len(self.Sy), len(self.Sz)
        self.N = [self.Nx, self.Ny, self.Nz]
        self.dtau = self.T / (self.Ntau)
        self.step = (Ntau * np.arange(int(self.T) + 1)).astype('int')   # [180, 360, 540, ...]
        
        self.u = np.zeros((self.Nx, self.Ny, self.Nz,))
        self.ki_u = np.zeros((self.Nx, self.Ny, self.Nz,))
        
    def _setup_coefficients(self):
        pass
    
    def _setup_boundary_conditions(self):
        pass
    
    def _traverse_time(self):
        pass
    
    def _redemption(self):
        pass
    
    def _traverse_grid(self):
        pass
    
    def _interpolate(self):
        pass
    
    def _get_grid(self):
        
        
        self._setup_coefficients()
        self._setup_boundary_conditions()
        self._traverse_time()
        
        
        
        return self.u, self.ki_u
    


class ELS3D(FDM3D):
    
    def _setup_coefficients(self):
        self.c = np.zeros((self.Nx, self.Ny, self.Nz))
        self.i_plus_one = np.zeros((self.Nx, self.Ny, self.Nz))
        self.j_plus_one = np.zeros((self.Nx, self.Ny, self.Nz))
        self.k_plus_one = np.zeros((self.Nx, self.Ny, self.Nz))
        self.i_minus_one = np.zeros((self.Nx, self.Ny, self.Nz))
        self.j_minus_one = np.zeros((self.Nx, self.Ny, self.Nz))
        self.k_minus_one = np.zeros((self.Nx, self.Ny, self.Nz))
        self.xy = np.zeros((self.Nx, self.Ny, self.Nz))
        self.yz = np.zeros((self.Nx, self.Ny, self.Nz))
        self.zx = np.zeros((self.Nx, self.Ny, self.Nz))
        self.i_j_k = np.zeros((self.Nx, self.Ny, self.Nz))

        for i in range(1, self.Nx - 1):
            for j in range(1, self.Ny - 1):
                for k in range(1 ,self.Nz - 1):
                    self.c[i][j][k] = 1 / self.dtau + 0.5 * self.r + (self.sigma_x**2 * self.Sx[i]**2 - self.r * self.Sx[i] * self.dSx[i]) / (self.dSx[i]*(self.dSx[i-1] + self.dSx[i])) +                                                                      (self.sigma_y**2 * self.Sy[j]**2 - self.r * self.Sy[j] * self.dSy[j]) / (self.dSy[j]*(self.dSy[j-1] + self.dSy[j])) +                                                                      (self.sigma_z**2 * self.Sz[k]**2 - self.r * self.Sz[k] * self.dSz[k]) / (self.dSz[k]*(self.dSz[k-1] + self.dSz[k]))
                    self.i_plus_one[i][j][k] = (self.sigma_x**2 * self.Sx[i]**2 + self.r * self.Sx[i] * self.dSx[i-1]) / (self.dSx[i]*(self.dSx[i-1] + self.dSx[i])) / self.c[i][j][k]
                    self.j_plus_one[i][j][k] = (self.sigma_y**2 * self.Sy[j]**2 + self.r * self.Sy[j] * self.dSy[j-1]) / (self.dSy[j]*(self.dSy[j-1] + self.dSy[j])) / self.c[i][j][k]
                    self.k_plus_one[i][j][k] = (self.sigma_z**2 * self.Sz[k]**2 + self.r * self.Sz[k] * self.dSz[k-1]) / (self.dSz[k]*(self.dSz[k-1] + self.dSz[k])) / self.c[i][j][k]
                    
                    self.i_minus_one[i][j][k] = (self.sigma_x**2 * self.Sx[i]**2 - self.r * self.Sx[i] * self.dSx[i]) / (self.dSx[i-1]*(self.dSx[i-1] + self.dSx[i])) / self.c[i][j][k]
                    self.j_minus_one[i][j][k] = (self.sigma_y**2 * self.Sy[j]**2 - self.r * self.Sy[j] * self.dSy[j]) / (self.dSy[j-1]*(self.dSy[j-1] + self.dSy[j])) / self.c[i][j][k]
                    self.k_minus_one[i][j][k] = (self.sigma_z**2 * self.Sz[k]**2 - self.r * self.Sz[k] * self.dSz[k]) / (self.dSz[k-1]*(self.dSz[k-1] + self.dSz[k])) / self.c[i][j][k]
                    
                    self.xy[i][j][k] = self.corr_xy * self.sigma_x * self.sigma_y * self.Sx[i] * self.Sy[j]                                         / (self.dSx[i] * self.dSy[j] + self.dSx[i-1] * self.dSy[j] + self.dSx[i] * self.dSy[j-1] + self.dSx[i-1] * self.dSy[j-1]) / self.c[i][j][k]
                    self.yz[i][j][k] = self.corr_yz * self.sigma_y * self.sigma_z * self.Sy[j] * self.Sz[k]                                         / (self.dSy[j] * self.dSz[k] + self.dSy[j-1] * self.dSz[k] + self.dSy[j] * self.dSz[k-1] + self.dSy[j-1] * self.dSz[k-1]) / self.c[i][j][k]
                    self.zx[i][j][k] = self.corr_zx * self.sigma_z * self.sigma_x * self.Sz[k] * self.Sx[i]                                         / (self.dSz[k] * self.dSx[i] + self.dSz[k-1] * self.dSx[i] + self.dSz[k] * self.dSx[i-1] + self.dSz[k-1] * self.dSx[i-1]) / self.c[i][j][k]
                    
                    self.i_j_k[i][j][k] = (1 / self.dtau - 0.5 * self.r) / self.c[i][j][k] - self.i_plus_one[i][j][k] - self.j_plus_one[i][j][k] - self.k_plus_one[i][j][k]
    
    
    def _setup_boundary_conditions(self):
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    if (self.Sx[i] <= self.barrier) or (self.Sx[j] <= self.barrier) or (self.Sx[k] <= self.barrier):
                        self.u[i, j, k] = self.face * np.min([self.Sx[i] / self.S0x, self.Sy[j] / self.S0y, self.Sz[k] / self.S0z])
                        self.ki_u[i, j, k] = self.face * np.min([self.Sx[i] / self.S0x, self.Sy[j] / self.S0y, self.Sz[k] / self.S0z])
                        
                    elif (self.Sx[i] < self.Kx[0]) or (self.Sy[j] < self.Ky[0]) or (self.Sz[k] < self.Kz[0]):
                        self.u[i, j, k] = self.face * (1 + self.dummy)
                        self.ki_u[i, j, k] = self.face * np.min([self.Sx[i] / self.S0x, self.Sy[j] / self.S0y, self.Sz[k] / self.S0z])
                    else:
                        self.u[i, j, k] = self.face * (1 + self.q[0])
                        self.ki_u[i, j, k] = self.face * (1 + self.q[0])


    def _traverse_time(self):
        n = 1
        idx_step = 1
        for tau in range(self.Ntau):
            if tau == self.step[idx_step]:
                self._redemption(idx_step)
                idx_step += 1
                #_printProgress(n, self.Ntau, 'Progress:', 'Complete', 0, 50)
                n += 1
            else:
                start = timeit.default_timer()
                self._traverse_grid()
                #_printProgress(n, self.Ntau, 'Progress:', 'Complete', 0, 50)
                stop = timeit.default_timer()
                time_delta = stop - start
                print('Runtime: ', time.strftime('%H:%M:%S', time.gmtime(time_delta)))
                n += 1
                

    def _redemption(self, idx_step):
        
        
        min_i = _get_argmin(self.Sx, self.Kx[idx_step])  
        min_j = _get_argmin(self.Sy, self.Ky[idx_step])  
        min_k = _get_argmin(self.Sz, self.Kz[idx_step])  
        
        self.u[min_i:, min_j:, min_k:] = self.face * (1 + self.q[idx_step])
        self.ki_u[min_i:, min_j:, min_k:] = self.face * (1 + self.q[idx_step])
        
        

    def _update_grid(self, i, j, k):
        self.u[i][j][k] = self.u[i+1][j][k] * self.i_plus_one[i][j][k] + self.u[i][j+1][k] * self.j_plus_one[i][j][k] + self.u[i][j][k+1] * self.k_plus_one[i][j][k] +                           self.u[i-1][j][k] * self.i_minus_one[i][j][k] + self.u[i][j-1][k] * self.j_minus_one[i][j][k] + self.u[i][j][k-1] * self.k_minus_one[i][j][k] +                           (self.u[i+1][j+1][k] - self.u[i-1][j+1][k] - self.u[i+1][j-1][k] + self.u[i-1][j-1][k]) * self.xy[i][j][k] +                           (self.u[i][j+1][k+1] - self.u[i][j-1][k+1] - self.u[i][j+1][k-1] + self.u[i][j-1][k-1]) * self.yz[i][j][k] +                           (self.u[i+1][j][k+1] - self.u[i+1][j][k-1] - self.u[i-1][j][k+1] + self.u[i-1][j][k-1]) * self.zx[i][j][k] +                           self.u[i][j][k] * self.i_j_k[i][j][k]
    
        self.ki_u[i][j][k] = self.ki_u[i+1][j][k] * self.i_plus_one[i][j][k] + self.ki_u[i][j+1][k] * self.j_plus_one[i][j][k] + self.ki_u[i][j][k+1] * self.k_plus_one[i][j][k] +                              self.ki_u[i-1][j][k] * self.i_minus_one[i][j][k] + self.ki_u[i][j-1][k] * self.j_minus_one[i][j][k] + self.ki_u[i][j][k-1] * self.k_minus_one[i][j][k] +                              (self.ki_u[i+1][j+1][k] - self.ki_u[i-1][j+1][k] - self.ki_u[i+1][j-1][k] + self.ki_u[i-1][j-1][k]) * self.xy[i][j][k] +                              (self.ki_u[i][j+1][k+1] - self.ki_u[i][j-1][k+1] - self.ki_u[i][j+1][k-1] + self.ki_u[i][j-1][k-1]) * self.yz[i][j][k] +                              (self.ki_u[i+1][j][k+1] - self.ki_u[i+1][j][k-1] - self.ki_u[i-1][j][k+1] + self.ki_u[i-1][j][k-1]) * self.zx[i][j][k] +                              self.ki_u[i][j][k] * self.i_j_k[i][j][k]
        
    def _neumann_boundary(self):
        self.u[self.Nx - 1, :self.Ny - 1, :self.Nz - 1] = self.u[self.Nx - 2, :self.Ny - 1, :self.Nz - 1]
        self.u[:self.Nx - 1, self.Ny - 1, :self.Nz - 1] = self.u[:self.Nx - 1, self.Ny - 2, :self.Nz - 1]
        self.u[:self.Nx - 1, :self.Ny - 1, self.Nz - 1] = self.u[:self.Nx - 1, :self.Ny - 1, self.Nz - 2]
        
        self.ki_u[self.Nx - 1, :self.Ny - 1, :self.Nz - 1] = self.ki_u[self.Nx - 2, :self.Ny - 1, :self.Nz - 1]
        self.ki_u[:self.Nx - 1, self.Ny - 1, :self.Nz - 1] = self.ki_u[:self.Nx - 1, self.Ny - 2, :self.Nz - 1]
        self.ki_u[:self.Nx - 1, :self.Ny - 1, self.Nz - 1] = self.ki_u[:self.Nx - 1, :self.Ny - 1, self.Nz - 2]
        

    def _traverse_grid(self):
        
        min_i = _get_argmin(self.Sx, self.barrier)
        min_j = _get_argmin(self.Sy, self.barrier)
        min_k = _get_argmin(self.Sz, self.barrier)
        
        self.u[:min_i, :, :] = copy.deepcopy(self.ki_u[:min_i :, :])
        self.u[:, :min_j, :] = copy.deepcopy(self.ki_u[:, :min_j, :])
        self.u[:, :, :min_k] = copy.deepcopy(self.ki_u[:, :, :min_k])
        
        
        for i in range(1, self.Nx - 1):
            for j in range(1, self.Ny - 1):
                for k in range(1, self.Nz - 1):
                    
                    self._update_grid(i, j, k)
                    
        self._neumann_boundary()
        
        for j in range(1, self.Ny - 1):
            for k in range(1, self.Nz - 1):
                for i in range(1, self.Nx - 1):
                    
                    self._update_grid(i, j, k)
                    
        self._neumann_boundary()
        
        for k in range(1, self.Nz - 1):
            for i in range(1, self.Nx - 1):
                for j in range(1, self.Ny - 1):
                    
                    self._update_grid(i, j, k)
                    
        self._neumann_boundary()
        
        for i in range(1, self.Nx - 1):
            for k in range(1, self.Nz - 1):
                for j in range(1, self.Ny - 1):
                    
                    self._update_grid(i, j, k)
                    
        self._neumann_boundary()
                    
        for j in range(1, self.Ny - 1):
            for i in range(1, self.Nx - 1):
                for k in range(1, self.Nz - 1):
                    
                    self._update_grid(i, j, k)
                    
        self._neumann_boundary()
        
        for k in range(1, self.Nz - 1):
            for j in range(1, self.Ny - 1):
                for i in range(1, self.Nx - 1):
                    
                    self._update_grid(i, j, k)

        self._neumann_boundary()
        
        
        
        
if __name__ == "__main___":
    els3d = ELS3D(S0 = [100, 100, 100], K = ([85, 85, 90, 90, 95, 95], [85, 85, 90, 90, 95, 95],[85, 85, 90, 90, 95, 95]),
                r = 0.03, T = 3, sigma = [0.3, 0.3, 0.3], Smax = [150, 150, 150], Ntau = 4320, dS = [15, 5, 20], corr = [0.5, 0.5, 0.5], 
                face = 100, q = [0.3, 0.25, 0.2, 0.15, 0.1, 0.05], dummy = 0.26, barrier = 50)

    grids = els3d._get_grid()


# In[ ]:





# In[21]:


els3d = ELS3D(S0 = [100, 100, 100], K = ([85, 85, 90, 90, 95, 95], [85, 85, 90, 90, 95, 95],[85, 85, 90, 90, 95, 95]),
                r = 0.03, T = 3, sigma = [0.3, 0.3, 0.3], Smax = [150, 150, 150], Ntau = 1540, dS = [5, 1, 5], corr = [0.5, 0.5, 0.5], 
                face = 100, q = [0.3, 0.25, 0.2, 0.15, 0.1, 0.05], dummy = 0.26, barrier = 50)

grids = els3d._get_grid()


# In[ ]:




