#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np


# In[8]:


class FDM:
    
    def __init__(self, S0, K, r, T, sigma, Smax, M, N, is_call=True):
        self.S0, self.K, self.r, self.T = S0, K, r, T
        self.sigma, self.Smax = sigma, Smax
        self.M, self.N = int(M), int(N)
        self.is_call = is_call
        
        self.dS = Smax / self.M
        self.dt = T / self.N
        self.i_values = np.arange(self.M)
        self.j_values = np.arange(self.N)
        self.grid = np.zeros((self.M + 1, self.N + 1))
        self.boundary_conds = np.linspace(0, self.Smax, self.M + 1)
        
    def _setup_boundary_conditions_(self):
        pass
    
    def _setup_coefficients_(self):
        pass
    
    def _traverse_grid_(self):
        # Iterate the grid backwards in time
        pass
    
    def _interpolate_(self):
        """
        Use piecewise linear interpolation on the initial
        grid column to get the closest price at S0
        """
        return np.interp(self.S0, self.boundary_conds, self.grid[:,0])
    
    def price(self):
        self._setup_boundary_conditions_()
        self._setup_coefficients_()
        self._traverse_grid_()
        return self._interpolate_()


# In[9]:


class FDMExplicitEU(FDM):
    
    def _setup_boundary_conditions_(self):
        if self.is_call:
            self.grid[:, -1] = np.maximum(self.boundary_conds - self.K, 0)
            self.grid[-1, :-1] = (self.Smax - self.K) * np.exp(-self.r * self.dt * (self.N - self.j_values))
            
        else:
            self.grid[:, -1] = np.maximum(self.K - self.boundary_conds, 0)
            self.grid[0, :-1] = (self.K) * np.exp(-self.r * self.dt * (self.N - self.j_values))
        
    def _setup_coefficients_(self):
        self.a = 0.5 * self.dt * (self.sigma**2 * self.i_values**2 - self.r * self.i_values)
        self.b = 1 - self.dt * (self.sigma**2 * self.i_values**2 + self.r)
        self.c = 0.5 * self.dt * (self.sigma**2 * self.i_values**2 + self.r * self.i_values)
    
    def _traverse_grid_(self):
        for j in reversed(self.j_values):
            for i in range(1, self.M):
                self.grid[i, j] = self.a[i] * self.grid[i - 1, j + 1] +                                   self.b[i] * self.grid[i, j + 1] +                                   self.c[i] * self.grid[i + 1, j + 1]
        

                


# In[10]:


option = FDMExplicitEU(50, 50, 0.1, 5/12, 0.4, 100, 150, 1500, False)
                    # S0, K, r, T, sigma, Smax, M, N, is_call=True
print(option.price())


# In[15]:


option.M


# In[13]:


len(option.a)


# In[11]:


import matplotlib.pyplot as plt
plt.plot(option.grid[:, 0])
plt.show()


# In[44]:


class FDMImplicitEU(FDM):
    
    def _setup_boundary_conditions_(self):
        if self.is_call:
            self.grid[:, -1] = np.maximum(self.boundary_conds - self.K, 0)
            self.grid[-1, :-1] = (self.Smax - self.K) * np.exp(-self.r * self.dt * (self.N - self.j_values))
            
        else:
            self.grid[:, -1] = np.maximum(self.K - self.boundary_conds, 0)
            self.grid[0, :-1] = (self.K) * np.exp(-self.r * self.dt * (self.N - self.j_values))
        
    def _setup_coefficients_(self):
        self.a = 0.5 * self.dt * (self.r * self.i_values - self.sigma**2 * self.i_values**2)
        self.b = 1 + self.dt * (self.sigma**2 * self.i_values**2 + self.r)
        self.c = - 0.5 * self.dt * (self.r * self.i_values + self.sigma**2 * self.i_values**2)
        
    def thomas(self, a, b, c, d):
    
        '''
        TDMA solver, a b c d can be NumPy array type or Python list type.
        refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        '''
        a = np.array(a, dtype=float)
        b = np.array(b, dtype=float)
        c = np.array(c, dtype=float)
        d = np.array(d, dtype=float)

        n = len(d) # the number of equations

        c[0] = c[0] / b[0]
        d[0] = d[0] / b[0]

        for i in range(1, n):
            c[i] = c[i] / (b[i] - a[i]*c[i-1])
            d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1])

        x = np.zeros(n)
        x[n-1] = d[n-1]

        for i in range(n-2, -1, -1):
            x[i] = d[i] - c[i] * x[i+1]

        return x
    
    def _traverse_grid_(self):
        for j in reversed(self.j_values):
            d = np.copy(self.grid[:, j + 1])
            d[1] -= self.a[1] * self.grid[0, j]
            d[-2] -= self.c[-2] * self.grid[-1, j]
            
            self.grid[1:-1, j] = self.thomas(self.a[1:], self.b[1:], self.c[1:], d[1:-1])

        


# In[45]:


option = FDMImplicitEU(50, 50, 0.1, 5/12, 0.4, 100, 100, 100, False)
option.price()


# In[ ]:




