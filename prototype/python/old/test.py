def get_psi(lat,lon):
  res = U_0*np.cos(lon*0.5)**2
  res = res*np.cos(lat)**2
  return res

delta = np.cos(lon_cur*0.5)*np.cos(lat_cur)/np.cos(lat_prev)
lon_prev = 2*np.arccos(delta)
 
