// functionTiming.cpp — MCMC progress reporting.
//
// Original code by Wayne Taylor from the bayesm package (Peter Rossi).

#include "maxdiff.h"

time_t itime;
char buf[100];

void startMcmcTimer() {
  itime = time(NULL);
  Rcout << " MCMC Iteration (est time to end - min) \n";
}

void infoMcmcTimer(int rep, int R) {
  time_t ctime = time(NULL);
  char buf[32];
  double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep + 1);
  snprintf(buf, 32, " %d (%.1f)\n", rep + 1, timetoend);
  Rcout << buf;
}

void endMcmcTimer() {
  time_t ctime = time(NULL);
  char buf[32];
  snprintf(buf, 32, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0);
  Rcout << buf;
  itime = 0;
}
