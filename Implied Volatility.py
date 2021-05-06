#!/usr/bin/env python
# coding: utf-8

# In[4]:


# Historical Volatility

import numpy as np

stock_price = np.array([230, 241, 238, 227, 243, 244, 253, 237, 239, 240, 253])
log_return = np.log(stock_price[1:] / stock_price[:-1])
mean_log_return = np.mean(log_return)
std_log_return = np.std(log_return, ddof=1)
annualied_volatility = np.sqrt(365) * std_log_return
annualied_volatility


# In[13]:


import numpy as np
from scipy.stats import norm

def get_bs_call(S, K, T, r, vol):
    d1 = (np.log(S / K) + T * (r + 0.5 * vol**2)) / (vol * np.sqrt(T))
    d2 = d1 - vol * np.sqrt(T)
    return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

def get_vega(S, K, T, r, vol):
    d1 = (np.log(S / K) + T * (r + 0.5 * vol**2)) / (vol * np.sqrt(T))
    return S * np.sqrt(T) * norm.pdf(d1)


# In[14]:


# Bisection Method

S = 100 ; K = 100 ; T = 1 ; C = 20 ; r = 0.05

vol_a = 0.0001 ; vol_b = 5
vol_mid = (vol_a + vol_b) / 2

tol = 1e-6
n = 0
while abs(vol_mid - vol_b) > tol:
    call_a = get_bs_call(S, K, T, r, vol_a)
    call_mid = get_bs_call(S, K, T, r, vol_mid)
    
    func_a = call_a - C
    func_mid = call_mid - C
    
    if func_a * func_mid <= 0 :
        vol_b = vol_mid
    else:
        vol_a = vol_mid
    
    vol_mid = (vol_a + vol_b) / 2
    n += 1

bisection_vol = vol_mid


# In[15]:


bisection_vol


# In[22]:


# Newton-Rapson Method

S = 100 ; K = 100 ; r = 0.05 ; T = 1 ; C = 20

vol = 0.2
tol = 1e-6
call = get_bs_call(S, K, T, r, vol)
vega = get_vega(S, K, T, r, vol)

while abs(call - C) > tol:
    vol -= (call - C) / vega
    
    call = get_bs_call(S, K, T, r, vol)
    vega = get_vega(S, K, T, r, vol)

implied_vol = vol


# In[23]:


implied_vol


# In[ ]:




