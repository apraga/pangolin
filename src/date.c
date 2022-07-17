/*-------------------------------------------------------------------------------
* CERFACS, Aviation and Environment
* -------------------------------------------------------------------------------
* 
*  MODULE: Date
* 
* > @author
* >  Alexis Praga
* 
*  DESCRIPTION: 
* 
* This converts a 32-bits integer date on the form YYYMMDDHHMM to a C date type
* It allows for date operation and the model can keep track of the date
* associated to the current iteration.
*
* Warning : this may not work for dates prior to 1900 (see C documentation).
*
* Note : for Fortran to call some functions, we need to use pointers only as
* arguments and add an underscore in the function name.
* */

#include <stdio.h>
#include <time.h>
int x = 0;

/* Start and end time */
time_t t_start, t_end;
struct tm start, end;
time_t dt;

void extract_date_(int *date, int *time, int *y, int *m, int *d, int *h, int *mi) {
  int dtmp = *date;
  int ttmp = *time;

  *mi = ttmp % 100;
  ttmp /= 100;
  *h = ttmp;

  *d = dtmp % 100;
  dtmp /= 100;
  *m = dtmp % 100;
  dtmp /= 100;
  *y = dtmp;
}


void set_time(int date, int time, struct tm *t, time_t *t_f) {
  t->tm_sec = 0;
  extract_date_(&date, &time, &t->tm_year, &t->tm_mon, &t->tm_mday, 
      &t->tm_hour, &t->tm_min);
  /* In respect to 1900 */
  t->tm_year = t->tm_year -1900;
  /* Normalize time (fill missing info)*/
  *t_f = mktime(t);
}

/* Input in seconds */
void set_timestep_date_(double *dt_loc) {
  dt = (*dt_loc);
}

/* Input as a date (YYYYMMDD) and time (HHMM) */
void set_startend_times_(int *d1, int *t1, int *d2, int *t2) {
  set_time(*d1, *t1, &start, &t_start);
  set_time(*d2, *t2, &end, &t_end);
  //printf("start time %d %d %lld \n", *d1, *t1, t_start);
}

/* dt in seconds */
void get_nb_iterations_(int* nb, double *dt) {
  (*nb) = (t_end - t_start)/(*dt);
}


/* Check if the time at the end of the iteration is a multiple of a period 
* Period is minutes */
void remainder_end_t_(int *r, int *iter, int *period) {
  int T = (int) (*period)*60;
  (*r) = ((*iter)*dt)  % T;
}

/* Get the date (YYYMMDD) and time (HHMM) at the beginning of a given iteration
 * starting from 1 */
void get_time_(int *date, int *time, int *iter) {
  time_t t_current = t_start + ((*iter)-1)*dt;
  struct tm *current = localtime(&t_current);
  int coef = 1;

  (*time) = current->tm_min + current->tm_hour*100;

  coef = 1;
  (*date) = current->tm_mday*coef;
  coef *= 100;
  (*date) += current->tm_mon*coef;
  coef *= 100;
  (*date) += (current->tm_year + 1900)*coef;
}

void print_t_start_() {
  printf("Start time %s \n",ctime (&t_start));
}
