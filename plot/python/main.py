# Script for calling any of the other python function
import error
import filament
import numpy as np
import datetime as d

#"filament_diag", "error_evol", "cfl_impact", "cv rate"
case = "error_evol"#filament_diag"
is_error = True
spec = "un" #"hourdin_norot"

# Create filament diagnostics
if case == "filament_diag":
  folder="/wkdir/pae2/praga/parallel/filaments/"
  tau = np.arange(0.1,1.1,0.05)
  lats = ['320']
  
  for lat in lats:
    print "lat2=", lat
    file0 = folder+"ratio_"+lat+"lat_"+spec+"_CFL1.0_0.dat"
    fileT2 = folder+"ratio_"+lat+"lat_"+spec+"_CFL1.0_T2.dat"
  
    l_f = filament.diagnostics(file0, fileT2, '1.5',tau)
    f = open("lf-pangolin-"+lat+"lat-"+spec+"-CFL1.0.dat", "w")
    for i in range(len(tau)):
      s = "%f %f\n" % (tau[i], l_f[i])
      f.write(s)
    f.close()

# Create error evolution
elif case == "error_evol":
  folder="/wkdir/pae2/praga/parallel/error/hourdin_80lat_1tracer/"
  
  l2 = []
  loo = []
  dates = []
  t = d.datetime(2013,01,03,00,00)
  dt = d.timedelta(hours=1)

  while t <= d.datetime(2013,01,03,00,00):
    s = "201301%02d%02d00" % (t.day,t.hour)
    print "t", t
    #file_num = folder+"ratio_1_"+s+".dat"
    file_num = "../../output/ratio_1_201301030000.dat"
    file_anal = folder+"ratio_anal_"+s+".dat"
    print file_anal
    l2.append(error.get_error(file_anal, file_num, '2'))
    loo.append(error.get_error(file_anal, file_num, 'oo'))
    dates.append(t.strftime("%Y-%m-%d %H:%M"))
    t += dt

  print "l2 loo", l2, loo
  f = open("error_80lat_1tracer.dat", "w")
  f.write("#date l2 loo\n")
  for i in range(len(loo)):
    s = "%s %f %f\n" % (dates[i], l2[i], loo[i])
    f.write(s)
  f.close()

# Computes errors for cfl impact
elif case == "cfl_impact":
  #dts = [0.8]#[24,20,19.2,18,16,12,10,8]
  dts = [9,7.5,6,5,4,3,0.6];

  l2 = []
  loo = []
  lquad = []
  folder="/wkdir/pae2/praga/parallel/cfl_impact/hourdin_80lat/"
  for dt in dts:
    file_num = folder+"ratio_dt_"+str(dt)+"_1_201301050000.dat"
    file_anal = folder+"ratio_anal_1_201301050000.dat"
    l2.append(error.get_error(file_anal, file_num, '2'))
    loo.append(error.get_error(file_anal, file_num, 'oo'))
    lquad.append(error.get_error(file_anal, file_num, 'quad'))

  print l2
  f = open("errors_norot.dat", "w")
  f.write("#dt l2 loo lquad\n")
  for i in range(len(loo)):
    s = "%f %f %f %f\n" % (dts[i], l2[i], loo[i], lquad[i])
    f.write(s)
  f.close()


  print loo

# Create convergence rates
elif case == "cv_rate":
  folder="/wkdir/pae2/praga/parallel/cv_rate/"
  
  l2 = []
  loo = []
  lats = [20, 40, 80, 160, 320]#, 640]
  for i in lats:
    file0 = folder+"ratio_"+str(i)+"lat_"+spec+"_CFL0.7_0.dat"
    fileT = folder+"ratio_"+str(i)+"lat_"+spec+"_CFL0.7_T.dat"
  
    l2.append(error.get_error(file0, fileT, '2'))
    loo.append(error.get_error(file0, fileT, 'oo'))
    
  f = open("conv-pangolin-"+spec+"-CFL0.7.dat", "w")
  res = [3, 1.5, 0.75, 0.375, 0.1875, 0.093]
  f.write("#resolution nb_cells l2 loo\n")
  for i in range(len(l2)):
    nb = 6*lats[i]**2
    s = "%f %d %f %f\n" % (res[i], nb, l2[i], loo[i])
    f.write(s)
  f.close()

# Create convergence rates for pole impact
elif case == "cv_rate_pole":
  folder="/wkdir/pae2/praga/parallel/cv_rate/"
  
  l2 = []
  loo = []
  lats = [20, 40, 80, 160, 320]#, 640]
  for i in lats:
    print i
    file0 = folder+"ratio_"+str(i)+"lat_"+spec+"_num_CFL0.7_T.dat"
    fileT = folder+"ratio_"+str(i)+"lat_"+spec+"_anal_T.dat"
  
    l2.append(error.get_error(file0, fileT, '2'))
    loo.append(error.get_error(file0, fileT, 'oo'))
  
  f = open("conv-pangolin-"+spec+"-CFL0.7.dat", "w")
  res = [3, 1.5, 0.75, 0.375, 0.1875, 0.093]
  f.write("#resolution nb_cells l2 loo\n")
  for i in range(len(l2)):
    nb = 6*lats[i]**2
    s = "%f %d %f %f\n" % (res[i], nb, l2[i], loo[i])
    f.write(s)
  f.close()


else:
  print "I don't understand what you want"
