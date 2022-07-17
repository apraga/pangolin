from pangolingrid import *
import shutil
from datetime import *

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
# Date as a string, time in seconds
def generate_winds(date, t1=None, t2=None):
  g.generate_u(folder+"u_"+date, test, t1, t2)
  g.generate_v(folder+"v_"+date, test, t1, t2)

# Date as a string
def correct_winds(date):
  # Eventually, find nb lat from ratio file
  #g.init_from_file(folder+"ratio_"+str(t_start)+".dat")

  ext = date
  u_old, u_new, v_old, v_new = [folder+"u_"+ext, folder+"u_corr_"+ext,
      folder+"v_"+ext, folder+"v_corr_"+ext]
  g.correct_winds(u_old, u_new, v_old, v_new)

  # Replace old files
  if (g.io_type != "hdf5" and replace):
    print("Replacing old ascii files")
    shutil.move(u_new+g.ext, u_old+g.ext)
    shutil.move(v_new+g.ext, v_old+g.ext)


#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

#winds_format = "center"
winds_format = "borders"
g = PangolinGrid(winds_format)

replace = True
nb_lat2 = 80
folder = "./"
case = "gaussian_hills"

#gen_winds = (case == "hourdin" or case == "custom") 
gen_winds = True

t_ref = datetime(2013,1,1,0,0)
t_start = t_ref
t_end = datetime(2013,1,1,0,0)

test = AnalyticalCase(case)#, True)
g.init(2*nb_lat2, "ascii")


date = datetime.strftime(t_start,"%Y%m%d%H%M")

# Generate time-dependent winds
time = timedelta(0)
dt = timedelta(days=2)#hours=1)

if (case == "hourdin"):
  test.set_t1(time.total_seconds())
date = datetime.strftime(t_start,"%Y%m%d%H%M")

t = t_start
while (t <= t_end):

  date = datetime.strftime(t,"%Y%m%d%H%M")
  g.generate_ratio(folder+"ratio_1_"+date, test)
  if (gen_winds):
    generate_winds(date, time.total_seconds(), (time+dt).total_seconds())
    correct_winds(date)
  t += dt
  time += dt

  if (case == "hourdin"):
    test.set_t2(time.total_seconds())

