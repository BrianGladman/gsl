#include <config.h>
#include <math.h>

#include "tests.h"

double book1 (double x) {
  return pow(x,alpha) * log(1/x) ;
}

double book2 (double x) {
  return pow(4.0,-alpha) / (pow((x-M_PI/4.0),2.0) + pow(16.0,-alpha)) ;
}

double book3 (double x) {
  return cos(pow(2.0,alpha) * sin(x)) ;
}

double book4 (double x) {
  return pow(fabs(x - (1.0/3.0)),alpha) ;
}

double book5 (double x) {
  return pow(fabs(x - (M_PI/4.0)),alpha) ;
}

double book6 (double x) {
  return 1 / ((x + 1 + pow(2.0,-alpha)) * sqrt(1-x*x)) ;
}

double book7 (double x) {
  return pow(sin(x), alpha-1) ;
}

double book8 (double x) {
  return pow(log(1/x), alpha-1) ;
}

double book9 (double x) {
  return exp(20*(x-1)) * sin(pow(2.0,alpha) * x) ;
}

double book10 (double x) {
  return 1/sqrt(x*(1-x)) ;
}

double book11 (double x) {
  return exp(-pow(2.0,-alpha)*x)*cos(x)/sqrt(x) ;
}

double book12 (double x) {
  return x*x * exp(-pow(2.0,-alpha)*x) ;
}

double book13 (double x) {
  return pow(x,alpha-1)/pow((1+10*x),2.0) ;
}

double book14 (double x) {
  return pow(2.0,-alpha)/(((x-1)*(x-1)+pow(4.0,-alpha))*(x-2)) ;
}


