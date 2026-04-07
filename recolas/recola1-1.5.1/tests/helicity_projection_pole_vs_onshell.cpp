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

  double p_full_prod[4][4], p_full_dec1[3][4], p_full_dec2[3][4];
  double A2[2], B2[1], C2[1], D2[1], E2[1];
  double M2prod, M2dec1, M2dec2, M2full;
  double MW, MZ, wwidth, alpha;
  int i;

  // Recola::set_output_file_rcl("*");
  // Recola::set_print_level_squared_amplitude_rcl(3);
  MW = 80.3579736098775;
  MZ = 91.1534806191828;
  wwidth = 2.08429899827822;
  alpha = 7.555310522369e-3;

  Recola::use_gfermi_scheme_and_set_alpha_rcl(alpha);
  Recola::set_pole_mass_w_rcl(MW,wwidth);
  Recola::set_pole_mass_z_rcl(MZ,0e0);
  Recola::set_complex_mass_scheme_rcl();

  Recola::set_resonant_particle_rcl("W-");
  Recola::define_process_rcl(1,"A A -> W-[-](nu_e~ e-) W+[-](mu+ nu_mu)","NLO");

  Recola::generate_processes_rcl();

  double p_full[6][4] = 
    {{5000.0000000000000, 0.0000000000000000, 0.0000000000000000, 5000.0000000000000},
     {5000.0000000000000, 0.0000000000000000, 0.0000000000000000,-5000.0000000000000},
     {3243.8985189361665, 695.23069130847887,-3166.9342850974826,-100.29516884221320},
     {1756.1014810638335, 433.87511031166298,-1700.8566566424388,-52.263122770570874},
     {2249.2998992648813,-494.85571955991441, 2193.9693904301166, 31.083235152524409},
     {2750.7001007351187,-634.25008206022721, 2673.8215513098039, 121.47505646025964}};

  Recola::compute_process_rcl(1, p_full, "NLO", A2);
  Recola::get_squared_amplitude_rcl(1,0,"LO",M2full);
  check_value(M2full, 4.8388318579782629e-05, 12);
  

  Recola::reset_recola_rcl();

  Recola::use_gfermi_scheme_and_set_alpha_rcl(alpha);
  Recola::set_pole_mass_w_rcl(MW,0e0);
  Recola::set_pole_mass_z_rcl(MZ,0e0);
  Recola::set_on_shell_scheme_rcl();

  Recola::define_process_rcl(1,"A A -> W-[-] W+[-]","NLO");
  Recola::define_process_rcl(2," W-[-] -> nu_e~ e-","NLO");
  Recola::define_process_rcl(3," W+[-] -> mu+ nu_mu","NLO");

  Recola::generate_processes_rcl();

  for (i=0; i<= 3; i = i + 1)
   {p_full_prod[0][i] = p_full[0][i];
    p_full_prod[1][i] = p_full[1][i];
    p_full_prod[2][i] = p_full[2][i] + p_full[3][i];
    p_full_prod[3][i] = p_full[4][i] + p_full[5][i];
    p_full_dec1[0][i] = p_full_prod[2][i];
    p_full_dec1[1][i] = p_full[2][i];
    p_full_dec1[2][i] = p_full[3][i];
    p_full_dec2[0][i] = p_full_prod[3][i];
    p_full_dec2[1][i] = p_full[4][i];
    p_full_dec2[2][i] = p_full[5][i];}

  Recola::compute_process_rcl(1, p_full_prod, "NLO", B2);
  Recola::compute_process_rcl(2, p_full_dec1, "NLO", C2);
  Recola::compute_process_rcl(3, p_full_dec2, "NLO", D2);
  Recola::get_squared_amplitude_rcl(1,0,"LO",M2prod);
  Recola::get_squared_amplitude_rcl(2,0,"LO",M2dec1);
  Recola::get_squared_amplitude_rcl(3,0,"LO",M2dec2);
  M2full = M2prod*M2dec1*M2dec2/(MW*wwidth)/(MW*wwidth)/(MW*wwidth)/(MW*wwidth);
  check_value(M2full, 4.8388318579781714e-05, 10);

  Recola::reset_recola_rcl();

  return 0;
} 
