/*
 * Examples of gnuplot_i.c usage
 * Compilation: gcc -Wall -g example.c gnuplot_i.c -o example -lm
 *
 * TODO
 * - need an example for gnuplot_splot_grid(), gnuplot_splot_obj(), gnuplot_plot_obj_xy()
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gnuplot_i.h"

#define SECONDS 1
#define NPOINTS 50

int main(int argc, char *argv[]) {
  gnuplot_ctrl *h1, *h2, *h3;
  double x[NPOINTS], y[NPOINTS], z[NPOINTS];
  int i;

  /** Initialize the gnuplot handle */

  printf("*** Examples of gnuplot control through C ***\n");
  h1 = gnuplot_init();

  /** Dumb terminal: ASCII-plot */
  /** Simplest usage of gnuplot_cmd function with fewest dependencies */

  printf("\n*** dumb terminal\n");
  gnuplot_setterm(h1, "dumb", 150, 40);
  gnuplot_cmd(h1, "plot sin(x) w lines, cos(x) w lines");
  sleep(SECONDS);

  /** Equations */

  printf("\n*** various equations\n");
  gnuplot_setterm(h1, "wxt", 900, 400);
  gnuplot_resetplot(h1);
  printf("y = sin(x)\n");
  gnuplot_plot_equation(h1, "sin(x)", "sine points");  // points is the default linestyle
  sleep(SECONDS);
  printf("y = sin(x)*cos(2*x)\n");
  gnuplot_setstyle(h1, "lines");
  gnuplot_plot_equation(h1, "sin(x)*cos(2*x)", "sines lines");
  sleep(SECONDS);
  printf("y = 2*x-1\n");
  gnuplot_setstyle(h1, "dots");
  gnuplot_plot_equation(h1, "2*x-1", "slope dots");
  sleep(SECONDS);

  /** Styles */

  printf("\n*** various styles\n");
  gnuplot_resetplot(h1);
  printf("sin(x) in linespoints\n");
  gnuplot_setstyle(h1, "linespoints");
  gnuplot_plot_equation(h1, "sin(x)", "sine linespoints");
  sleep(SECONDS);
  printf("cos(x/pi) in impulses\n");
  gnuplot_setstyle(h1, "impulses");
  // note that M_PI cannot be used in the denominator, since that is not known to gnuplot...
  gnuplot_plot_equation(h1, "cos(x/(atan(1)*4))", "sine impulses");
  sleep(SECONDS);
  printf("atan(x/pi) in steps\n");
  gnuplot_setstyle(h1, "steps");
  // ...but it is possible to use gnuplot constants, such as "pi"
  gnuplot_plot_equation(h1, "atan(x/pi)", "arctangens steps");
  sleep(SECONDS);

  /** User defined 1d and 2d point sets */

  printf("\n*** user-defined list of points: 1d\n");
  gnuplot_resetplot(h1);
  printf("quadratic 1d array\n");
  gnuplot_setstyle(h1, "points");
  gnuplot_set_axislabel(h1, "x", "X");
  for (i = 0; i < NPOINTS; i++) {
    x[i] = (double)i * i;
  }
  gnuplot_plot_coordinates(h1, x, y, NPOINTS, "user-defined points: 1d");   // note: y is NULL
  sleep(SECONDS);

  printf("\n*** user-defined lists of points: 2d\n");
  gnuplot_resetplot(h1);
  printf("square root 2d array\n");
  gnuplot_setstyle(h1, "points");
  gnuplot_set_axislabel(h1, "y", "square root");
  for (i = 0; i < NPOINTS; i++) {
    x[i] = (double)i/2;  // change of axis scale
    y[i] = sqrt(i);
  }
  gnuplot_plot_coordinates(h1, x, y, NPOINTS, "user-defined points: 2d");
  sleep(SECONDS);

  /** gnuplot_plot_once demo: User defined 1d and 2d point sets */

  printf("\n*** gnuplot_plot_once: user-defined list of doubles\n");
  gnuplot_plot_once ("lines", "X", "Y", x, y, NPOINTS, "list of doubles (plot once)");

  /** Splot (surface plot) example */

  printf("\n*** parametric 3D plot\n");
  gnuplot_resetplot(h1);
  printf("lissajous curve in 3d space\n");
  gnuplot_setstyle(h1, "lines");
  gnuplot_set_axislabel(h1, "x", "X");
  gnuplot_set_axislabel(h1, "y", "Y");
  gnuplot_set_axislabel(h1, "z", "Z-axis");
  for (int i = 0; i < NPOINTS; i++) {
    x[i] = 2*sin((double)i/3);
    y[i] = 5*sin((double)i/2 + 1);
    z[i] = x[i] + y[i];
  }
  printf("Internal debugging information:\n");
  print_gnuplot_handle(h1);
  gnuplot_splot(h1, x, y, z, NPOINTS, "Lissajous");
  sleep(SECONDS);

  /** Contour plot example */

  printf("\n*** contour plot\n");
  gnuplot_resetplot(h1);
  printf("contour of multivariate bell curve\n");
  double xx[NPOINTS*NPOINTS], yy[NPOINTS*NPOINTS], zz[NPOINTS*NPOINTS];
  for (int i = 0; i < NPOINTS; i++) {
    for (int j = 0; j < NPOINTS; j++) {
      xx[NPOINTS*i+j] = i;
      yy[NPOINTS*i+j] = j;
      zz[NPOINTS*i+j] = 1000*sqrt(pow(i-NPOINTS/2, 2)+pow(j-NPOINTS/2, 2));
    }
  }
  gnuplot_setstyle(h1, "lines");
  gnuplot_contour_plot(h1, xx, yy, zz, NPOINTS, NPOINTS, "Contour");
  sleep(SECONDS);

  /** Scatter plot: gnuplot example with data file */
  /** Note that the program exits gracefully if file is not present, not readable or mispelled **/

  printf("\n*** scatter plot: data file\n");
  gnuplot_resetplot(h1);
  gnuplot_set_axislabel(h1, "x", "avg(mRS)");
  gnuplot_set_axislabel(h1, "y", "alpha");
  gnuplot_cmd(h1, "plot 'scatter.data' with points pointtype '•' linecolor 'blue'");
  sleep(SECONDS);

  /** Contour plot with color palette*/
  /** if no color palette is found, then gnuplots standard color palette will be used: rgbformulae */
  /** This does not work yet: no z-axis is shown */

  printf("\n*** contour plot with palette\n");
  gnuplot_resetplot(h1);
  gnuplot_set_axislabel(h1, "x", "shape");
  gnuplot_set_axislabel(h1, "y", "scale");
  gnuplot_set_axislabel(h1, "z", "alpha");
  gnuplot_cmd(h1, "set parametric");
  gnuplot_cmd(h1, "set contour base");
  gnuplot_cmd(h1, "set style data lines");
  gnuplot_cmd(h1, "set dgrid3d 20, 20, 20");
  gnuplot_cmd(h1, "set title 'Alpha-parameter by shape and scale'");
  gnuplot_cmd(h1, "load 'bentcoolwarm.palette'");
  gnuplot_cmd(h1, "set pm3d");
  gnuplot_cmd(h1, "splot 'indicatorvalues.data'");
  sleep(SECONDS);

  /** Multiple output screens */

  printf("\n*** multiple output windows\n");
  gnuplot_resetplot(h1);
  gnuplot_setstyle(h1, "lines");
  h2 = gnuplot_init();
  gnuplot_setstyle(h2, "lines");
  h3 = gnuplot_init();
  gnuplot_setstyle(h3, "filledcurves");
  printf("window 1: x*sin(x)\n");
  gnuplot_cmd(h1, "set linecolor rgb 'blue'");   // this does not work yet
  gnuplot_plot_equation(h1, "x*sin(x)", "x*sin(x)");
  sleep(SECONDS);
  printf("window 2: log(x)/x\n");
  gnuplot_plot_equation(h2, "log(x)/x", "log(x)/x");
  sleep(SECONDS);
  printf("window 3: sin(x)/x\n");
  gnuplot_plot_equation(h3, "sin(x)/x", "sin(x)/x");
  sleep(SECONDS);

  /** Close gnuplot handles */

  printf("\n*** closing all gnuplot windows\n");
  gnuplot_close(h1);
  gnuplot_close(h2);
  gnuplot_close(h3);
  return 0;
}
