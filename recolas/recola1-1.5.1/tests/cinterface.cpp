//#####################################################################
//
//  File   cinterface.cpp
//  is part of RECOLA (REcursive Computation of One Loop Amplitudes)
//
//  Copyright (C) 2015-2017   Stefano Actis, Ansgar Denner,
//                            Lars Hofer, Jean-Nicolas Lang,
//                            Andreas Scharf, Sandro Uccirati
//
//  RECOLA is licenced under the GNU GPL version 3,
//         see COPYING for details.
//
//#####################################################################
//
#include <complex>
#include "recola.h"
#include <cmath>
#include <iostream>

void check_value(double v1, double v2, int digits){
    double acc = pow(10.,-digits);
    double diff = fabs(v1-v2) / fabs(v1+v2) * 2;
    if (diff > acc) {
        std::cout << "Error: " << v1 << " != " << v2 << " (diff = " << diff << ")" << std::endl;
        exit(1);
    }
}

int main(int argc, char *argv[])
{

  Recola::set_print_level_squared_amplitude_rcl (1);

  Recola::set_print_level_correlations_rcl (1);

  Recola::define_process_rcl(1,"g g -> d d~ e+ e-","LO");
  Recola::define_process_rcl(2,"u g -> d d~ e+ e- u","LO");
  Recola::define_process_rcl(3,"A A -> mu+ mu- e+ e-","LO");
  Recola::define_process_rcl(4,"u A -> mu+ mu- e+ e- u","LO");

  Recola::unselect_all_gs_powers_BornAmpl_rcl(1);
  Recola::select_gs_power_BornAmpl_rcl(1,2);

  Recola::unselect_all_gs_powers_BornAmpl_rcl(2);
  Recola::select_gs_power_BornAmpl_rcl(2,3);

  Recola::unselect_all_gs_powers_BornAmpl_rcl(3);
  Recola::select_gs_power_BornAmpl_rcl(3,0);

  Recola::unselect_all_gs_powers_BornAmpl_rcl(4);
  Recola::select_gs_power_BornAmpl_rcl(4,0);

  Recola::generate_processes_rcl();

  double p[6][4] =
   {{2789.36556449,  0.00000000000,  0.00000000000,  2789.36556449},
    {2003.32704474,  0.00000000000,  0.00000000000, -2003.32704474},
    {1482.88702577,  87.3273179083, -1127.68979373,  958.980500246},
    {2032.83040303, -234.591543535,  1812.27464017,  890.520568996},
    {1220.45323502,  159.141025786, -661.270470691, -1013.36153340},
    {56.5219454083, -11.8768001589, -23.3143757501, -50.1010160985}};

  double k[7][4] =
   {{2914.11325883,  0.00000000000,  0.00000000000,  2914.11325883},
    {2003.32704474,  0.00000000000,  0.00000000000, -2003.32704474},
    {1321.46369109,  169.108121020, -876.002925898,  974.826961102},
    {2032.83040303, -234.591543535,  1812.27464017,  890.520568996},
    {1220.45323502,  159.141025786, -661.270470691, -1013.36153340},
    {56.5219454083, -11.8768001589, -23.3143757501, -50.1010160985},
    {286.171029020, -81.7808031113, -251.686867829,  108.901233489}};


  const int nlegs = Recola::get_legs_rcl(1);
  const int ncols = Recola::get_cs_rcl(1);
  int cmatrix[nlegs*ncols];

  Recola::get_colour_configurations_rcl (1,cmatrix);
  std::cout << "Colour configurations for process 1:" << std::endl;
  for (int i = 0; i < ncols; ++i) {
    for (int j = 0; j < nlegs; ++j) {
      std::cout << cmatrix[nlegs*i + j] << " ";
    }
    std::cout << std::endl;
  }

  const int nhels = Recola::get_cf_rcl(1);

  int hmatrix[nlegs*nhels];
  std::cout << "Helicity configurations for process 1:" << std::endl;
  Recola::get_helicity_configurations_rcl (1,hmatrix);
  for (int i = 0; i < nhels; ++i) {
    for (int j = 0; j < nlegs; ++j) {
      std::cout << hmatrix[nlegs*i + j] << " ";
    }
    std::cout << std::endl;
  }

  double A2;
  Recola::compute_process_rcl(1,p,"LO");
  Recola::get_squared_amplitude_rcl(1, 2, "LO", A2);
  check_value(A2, 0.26677515896114E-10, 10);

  Recola::compute_process_rcl(2,k,"LO");
  Recola::get_squared_amplitude_rcl(2, 3, "LO", A2);
  check_value(A2, 0.21229780305249E-15, 10);

  Recola::compute_process_rcl(3,p,"LO");
  Recola::get_squared_amplitude_rcl(3, 0, "LO", A2);
  check_value(A2, 0.30099594967041E-11, 10);

  Recola::compute_process_rcl(4,k,"LO");
  Recola::get_squared_amplitude_rcl(4, 0, "LO", A2);
  check_value(A2, 0.40093279151328E-18, 10);

  Recola::compute_colour_correlation_rcl(1,p,1,3);
  std::complex<double> v[4] =
    {std::complex<double>(-1151.5040504618346, 0.),
     std::complex<double>(-497.28644383654944, 0.),
     std::complex<double>( 580.15008908958009, 0.),
     std::complex<double>(-1151.5040504618348, 0.)};
  Recola::compute_spin_colour_correlation_rcl(1,p,1,3,v);
  Recola::get_spin_colour_correlation_rcl(1, 2, 1, 3, A2);
  check_value(A2, -0.16769090531220E-05, 10);

  Recola::compute_spin_correlation_rcl(3,p,1,v);
  Recola::get_spin_correlation_rcl(3, 0, A2);
  check_value(A2, 0.78379130944946E-06, 10);

  Recola::reset_recola_rcl();

  return 0;
}
//#####################################################################
