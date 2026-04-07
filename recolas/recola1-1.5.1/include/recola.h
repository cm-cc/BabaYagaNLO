//####################################################################
//
//  File  recola.h
//  is part of RECOLA (REcursive Computation of One Loop Amplitudes)
//
//  Copyright (C) 2015, 2016   Stefano Actis, Ansgar Denner,
//                             Lars Hofer, Jean-Nicolas Lang,
//                             Andreas Scharf, Sandro Uccirati
//
//  C++ interface developed by Benedikt Biedermann and Mathieu Pellen
//
//  RECOLA is licenced under the GNU GPL version 3,
//         see COPYING for details.
//
//####################################################################

// Macro to add name-mangling bits to fortran symbols.

// Define GNU_NAME_MANGLING(INTEL_NAME_MANGLING) if Recola has been build with
// gfortran(ifort)

#if defined(GNU_NAME_MANGLING)
// gfortran
//#pragma message ("Name mangling according to GNU." )
#define MODFUNCNAME(mod,fname) __ ## mod ## _MOD_ ## fname
#elif defined(INTEL_NAME_MANGLING)
// ifort
//#pragma message ("Name mangling according to INTEL." )
#define MODFUNCNAME(mod,fname) mod ## _mp_ ## fname ## _
#else
#warning "Name mangling not set. Trying GNU."
#define MODFUNCNAME(mod,fname) __ ## mod ## _MOD_ ## fname
#endif // if __GNUC__

#define RCLMOD(mod,fname) MODFUNCNAME(mod ## _rcl, fname ## _rcl)

#ifndef RECOLA_H_
#define RECOLA_H_

#include <complex>
struct dcomplex{double dr,di;};

#ifdef __cplusplus
extern "C" {
#endif

  // module input_rcl
  void RCLMOD(input,set_pole_mass_z)
       (const double*,const double*);
  void RCLMOD(input,set_onshell_mass_z)
       (const double*,const double*);
  void RCLMOD(input,set_pole_mass_w)
       (const double*,const double*);
  void RCLMOD(input,set_onshell_mass_w)
       (const double*,const double*);
  void RCLMOD(input,set_pole_mass_h)
       (const double*,const double*);
  void RCLMOD(input,set_pole_mass_electron)
       (const double*);
  void RCLMOD(input,set_pole_mass_muon)
       (const double*,const double*);
  void RCLMOD(input,set_pole_mass_tau)
       (const double*,const double*);
  void RCLMOD(input,set_pole_mass_up)
       (const double*);
  void RCLMOD(input,set_pole_mass_down)
       (const double*);
  void RCLMOD(input,set_pole_mass_charm)
       (const double*,const double*);
  void RCLMOD(input,set_pole_mass_strange)
       (const double*);
  void RCLMOD(input,set_pole_mass_top)
       (const double*,const double*);
  void RCLMOD(input,set_pole_mass_bottom)
       (const double*,const double*);
  void RCLMOD(input,set_light_fermions)
       (const double*);
  void RCLMOD(input,set_light_electron)
       ();
  void RCLMOD(input,set_light_muon)
       ();
  void RCLMOD(input,set_light_tau)
       ();
  void RCLMOD(input,set_light_up)
       ();
  void RCLMOD(input,set_light_down)
       ();
  void RCLMOD(input,set_light_charm)
       ();
  void RCLMOD(input,set_light_strange)
       ();
  void RCLMOD(input,set_light_top)
       ();
  void RCLMOD(input,set_light_bottom)
       ();
  void RCLMOD(input,unset_light_electron)
       ();
  void RCLMOD(input,unset_light_muon)
       ();
  void RCLMOD(input,unset_light_tau)
       ();
  void RCLMOD(input,unset_light_up)
       ();
  void RCLMOD(input,unset_light_down)
       ();
  void RCLMOD(input,unset_light_charm)
       ();
  void RCLMOD(input,unset_light_strange)
       ();
  void RCLMOD(input,unset_light_top)
       ();
  void RCLMOD(input,unset_light_bottom)
       ();
  void RCLMOD(input,use_dim_reg_soft)
       ();
  void RCLMOD(input,use_mass_reg_soft)
       (const double*);
  void RCLMOD(input,set_mass_reg_soft)
       (const double*);
  void RCLMOD(input,set_delta_uv)
       (const double*);
  void RCLMOD(input,set_mu_uv)
       (const double*);
  void RCLMOD(input,set_delta_ir)
       (const double*, const double*);
  void RCLMOD(input,set_mu_ir)
       (const double*);
  void RCLMOD(input,set_alphas)
       (const double*,const double*,const int*);
  void RCLMOD(input,get_alphas)
       (double*);
  void RCLMOD(input,get_renormalization_scale)
       (double*);
  void RCLMOD(input,get_flavour_scheme)
       (int*);
  void RCLMOD(input,set_alphas_masses)
       (const double*,const double*,const double*,
        const double*,const double*,const double*);
  void RCLMOD(wrapper,set_alphas_masses_nowidtharg)
       (const double*,const double*,const double*);
  void RCLMOD(wrapper,use_gfermi_scheme_noarg)
       ();
  void RCLMOD(wrapper,use_gfermi_scheme_and_set_alpha)
       (const double*);
  void RCLMOD(wrapper,use_gfermi_scheme_and_set_gfermi)
       (const double*);
  void RCLMOD(wrapper,use_gfermi_scheme_real_and_set_gfermi)
       (const double*);
  void RCLMOD(wrapper,use_gfermi_scheme_complex_and_set_gfermi)
       (const double*);
  void RCLMOD(input,use_alpha0_scheme)
       (const double*);
  void RCLMOD(wrapper,use_alpha0_scheme_noarg)
       ();
  void RCLMOD(input,use_alphaz_scheme)
       (const double*);
  void RCLMOD(wrapper,use_alphaz_scheme_noarg)
       ();
  void RCLMOD(input,get_alpha)
       (double*);
  void RCLMOD(input,set_complex_mass_scheme)
       ();
  void RCLMOD(input,set_on_shell_scheme)
       ();
  void RCLMOD(input,set_resonant_particle)
       (const char*,long int);
  void RCLMOD(input,set_quarkline)
       (const int*,const int*,const int*);
  void RCLMOD(input,switchon_resonant_selfenergies)
       ();
  void RCLMOD(input,switchoff_resonant_selfenergies)
       ();
  void RCLMOD(input,set_dynamic_settings)
       (const int*);
  void RCLMOD(input,set_momenta_correction)
       (const int*);
  void RCLMOD(input,set_draw_level_branches)
       (const int*);
  void RCLMOD(input,set_print_level_amplitude)
       (const int*);
  void RCLMOD(input,set_print_level_squared_amplitude)
       (const int*);
  void RCLMOD(input,set_print_level_correlations)
       (const int*);
  void RCLMOD(input,set_print_level_ram)
       (const int*);
  void RCLMOD(input,scale_coupling3)
       (dcomplex*,const char*,const char*, const char*,
        long int,long int,long int);
  void RCLMOD(input,scale_coupling4)
       (dcomplex*, const char*, const char*, const char*,const char*,
        long int,long int,long int,long int);
  void RCLMOD(input,switchoff_coupling3)
       (const char*,const char*,const char*,long int,long int,long int);
  void RCLMOD(input,switchoff_coupling4)
       (const char*,const char*,const char*,const char*,
        long int,long int,long int,long int);
  void RCLMOD(input,set_ifail)
       (const int*);
  void RCLMOD(input,get_ifail)
       (int*);
  void RCLMOD(input,set_collier_output_dir)
       (const char*,long int);
  void RCLMOD(input,set_compute_ir_poles)
       (const int*);
  void RCLMOD(input,set_output_file)
       (const char*,long int);

  // module process_definition_rcl
  void RCLMOD(process_definition,define_process)
       (const int*,const char*,const char*,long int,long int);
  void RCLMOD(wrapper,wrapper_set_gs_power)
       (const int*, const int[][2], const int*);
  void RCLMOD(process_definition,select_gs_power_bornampl)
       (const int*, const int*);
  void RCLMOD(process_definition,select_gs_power_loopampl)
       (const int*, const int*);
  void RCLMOD(process_definition,unselect_gs_power_bornampl)
       (const int*, const int*);
  void RCLMOD(process_definition,unselect_gs_power_loopampl)
       (const int*, const int*);
  void RCLMOD(process_definition,select_all_gs_powers_bornampl)
       (const int*);
  void RCLMOD(process_definition,select_all_gs_powers_loopampl)
       (const int*);
  void RCLMOD(process_definition,unselect_all_gs_powers_bornampl)
       (const int*);
  void RCLMOD(process_definition,unselect_all_gs_powers_loopampl)
       (const int*);
  void RCLMOD(process_definition,split_collier_cache)
       (const int*,const int*);

  // module process_generation_rcl
  void RCLMOD(process_generation,generate_processes)
       ();
  void RCLMOD(process_generation,process_exists)
       (const int*, int *);

  // module process_computation_rcl
  void RCLMOD(process_computation,set_resonant_squared_momentum)
       (const int*, const int*, const double*);
  void RCLMOD(process_computation,compute_running_alphas)
       (const double*, const int*, const int*);
  void RCLMOD(wrapper,wrapper_set_outgoing_momenta)
       (const int*, double[][4], double[][4], const int*);
  void RCLMOD(wrapper,wrapper_compute_process)
       (const int*, const double[][4], const int*,
        const char*, double*, int*, long int);
  void RCLMOD(wrapper,wrapper_rescale_process)
       (const int*, const char*, double*, long int);
  void RCLMOD(wrapper,wrapper_get_colour_configurations)
       (int*,int*,int*,int*);

  void RCLMOD(wrapper,wrapper_get_helicity_configurations)
       (int*,int*,int*,int*);
  // get_colour_configurations_rcl

  // get_helicity_configurations_rcl

  void RCLMOD(wrapper,wrapper_get_amplitude)
       (const int*, const int*, const char*, const int[], const int[],
        const int*, dcomplex*, long int);
  void RCLMOD(process_computation,get_squared_amplitude)
       (const int*, const int*, const char*, double*, long int);
  void RCLMOD(wrapper,wrapper_get_polarized_squared_amplitude)
       (const int*, const int*, const char*, const int[],
        const int*, double*, long int);
  void RCLMOD(wrapper,wrapper_compute_colour_correlation)
       (const int*, const double[][4], const int*, const int*,
        const int*, double*, int*);
  void RCLMOD(wrapper,wrapper_compute_colour_correlation_int)
       (const int*, const double[][4], const int*, const int*,
        const int*, double*, int*);
  void RCLMOD(wrapper,wrapper_compute_all_colour_correlations)
       (const int*, const double[][4], const int*, int*);
  void RCLMOD(wrapper,wrapper_compute_all_colour_correlations_int)
       (const int*, const double[][4], const int*, int*);
  void RCLMOD(wrapper,wrapper_rescale_colour_correlation)
       (const int*, const int*, const int*, double*);
  void RCLMOD(wrapper,wrapper_rescale_colour_correlation_int)
       (const int*, const int*, const int*, double*);
  void RCLMOD(process_computation,rescale_all_colour_correlations)
       (const int*);
  void RCLMOD(process_computation,get_colour_correlation)
       (const int*, const int*, const int*, const int*, double*);
  void RCLMOD(process_computation,get_colour_correlation_int)
       (const int*, const int*, const int*, const int*, double*);
  void RCLMOD(wrapper,wrapper_compute_spin_colour_correlation)
       (const int*, const double [][4], const int*, const int*, const int*,
        const dcomplex[4], double*, int*);
  void RCLMOD(wrapper,wrapper_rescale_spin_colour_correlation)
       (const int*, const int*, const int*,const dcomplex[4], double*);
  void RCLMOD(process_computation,get_spin_colour_correlation)
       (const int*, const int*, const int*, const int*, double*);
  void RCLMOD(wrapper,wrapper_compute_spin_correlation)
       (const int*, const double[][4], const int*, const int*,
        const dcomplex[4], double*, int*);
  void RCLMOD(wrapper,wrapper_rescale_spin_correlation)
       (const int*, const int*, const dcomplex[4], double*);
  void RCLMOD(process_computation,get_spin_correlation)
       (const int*, const int*, double*);
  void RCLMOD(wrapper,wrapper_compute_spin_correlation_matrix)
       (const int*, const double[][4], const int*, const int*,
        double[4][4], int*);
  void RCLMOD(wrapper,wrapper_rescale_spin_correlation_matrix)
       (const int*, const int*, double[4][4]);
  void RCLMOD(process_computation,get_spin_correlation_matrix)
       (const int*, const int*, double[4][4]);
  void RCLMOD(wrapper,get_legs)
       (const int*, int*);
  void RCLMOD(wrapper,get_cs)
       (const int*, int*);
  void RCLMOD(wrapper,get_cf)
       (const int*, int*);
  void RCLMOD(wrapper,wrapper_get_momenta)
       (const int*, double [][4], const int*);
  void RCLMOD(wrapper,wrapper_get_recola_version)
     (char*,int*,long int);
  void RCLMOD(process_computation,set_tis_required_accuracy)
       (const double*);
  void RCLMOD(process_computation,get_tis_required_accuracy)
       (double*);
  void RCLMOD(process_computation,set_tis_critical_accuracy)
       (const double*);
  void RCLMOD(process_computation,get_tis_critical_accuracy)
       (double*);
  void RCLMOD(process_computation,get_tis_accuracy_flag)
       (int*);

  // module reset_rcl
  void RCLMOD(reset,reset_recola)
       ();

#ifdef __cplusplus
}
#endif

namespace Recola{
/*
  module input_rcl
*/
inline void set_pole_mass_z_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_z)(&m,&g);
}
inline void set_onshell_mass_z_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_onshell_mass_z)(&m,&g);
}
inline void set_pole_mass_w_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_w)(&m,&g);
}
inline void set_onshell_mass_w_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_onshell_mass_w)(&m,&g);
}
inline void set_pole_mass_h_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_h)(&m,&g);
}
inline void set_pole_mass_electron_rcl
            (const double m)
{
  RCLMOD(input,set_pole_mass_electron)(&m);
}
inline void set_pole_mass_muon_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_muon)(&m,&g);
}
inline void set_pole_mass_tau_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_tau)(&m,&g);
}
inline void set_pole_mass_up_rcl
            (const double m)
{
  RCLMOD(input,set_pole_mass_up)(&m);
}
inline void set_pole_mass_down_rcl
            (const double m)
{
  RCLMOD(input,set_pole_mass_down)(&m);
}
inline void set_pole_mass_charm_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_charm)(&m,&g);
}
inline void set_pole_mass_strange_rcl
            (const double m)
{
  RCLMOD(input,set_pole_mass_strange)(&m);
}
inline void set_pole_mass_top_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_top)(&m,&g);
}
inline void set_pole_mass_bottom_rcl
            (const double m, const double g)
{
  RCLMOD(input,set_pole_mass_bottom)(&m,&g);
}
inline void set_light_fermions_rcl
            (const double m)
{
  RCLMOD(input,set_light_fermions)(&m);
}
inline void set_light_electron_rcl()
{
  RCLMOD(input,set_light_electron)();
}
inline void set_light_muon_rcl()
{
  RCLMOD(input,set_light_muon)();
}
inline void set_light_tau_rcl()
{
  RCLMOD(input,set_light_tau)();
}
inline void set_light_up_rcl()
{
  RCLMOD(input,set_light_up)();
}
inline void set_light_down_rcl()
{
  RCLMOD(input,set_light_down)();
}
inline void set_light_charm_rcl()
{
  RCLMOD(input,set_light_charm)();
}
inline void set_light_strange_rcl()
{
  RCLMOD(input,set_light_strange)();
}
inline void set_light_top_rcl()
{
  RCLMOD(input,set_light_top)();
}
inline void set_light_bottom_rcl()
{
  RCLMOD(input,set_light_bottom)();
}
inline void unset_light_electron_rcl()
{
  RCLMOD(input,unset_light_electron)();
}
inline void unset_light_muon_rcl()
{
  RCLMOD(input,unset_light_muon)();
}
inline void unset_light_tau_rcl()
{
  RCLMOD(input,unset_light_tau)();
}
inline void unset_light_up_rcl()
{
  RCLMOD(input,unset_light_up)();
}
inline void unset_light_down_rcl()
{
  RCLMOD(input,unset_light_down)();
}
inline void unset_light_charm_rcl()
{
  RCLMOD(input,unset_light_charm)();
}
inline void unset_light_strange_rcl()
{
  RCLMOD(input,unset_light_strange)();
}
inline void unset_light_top_rcl()
{
  RCLMOD(input,unset_light_top)();
}
inline void unset_light_bottom_rcl()
{
  RCLMOD(input,unset_light_bottom)();
}
inline void use_dim_reg_soft_rcl()
{
  RCLMOD(input,use_dim_reg_soft)();
}
inline void use_mass_reg_soft_rcl
            (const double m)
{
  RCLMOD(input,use_mass_reg_soft)(&m);
}
inline void set_mass_reg_soft_rcl
            (const double m)
{
  RCLMOD(input,set_mass_reg_soft)(&m);
}
inline void set_delta_uv_rcl
            (const double d)
{
  RCLMOD(input,set_delta_uv)(&d);
}
inline void set_mu_uv_rcl
            (const double m)
{
  RCLMOD(input,set_mu_uv)(&m);
}
inline void set_delta_ir_rcl
            (const double d1, const double d2)
{
  RCLMOD(input,set_delta_ir)(&d1,&d2);
}
inline void set_mu_ir_rcl
            (const double m)
{
  RCLMOD(input,set_mu_ir)(&m);
}
inline void set_alphas_rcl
            (const double a, const double s, const int nf)
{
  RCLMOD(input,set_alphas)(&a,&s,&nf);
}
inline void get_alphas_rcl
            (double &a)
{
  RCLMOD(input,get_alphas)(&a);
}
inline void get_renormalization_scale_rcl
            (double &mu)
{
  RCLMOD(input,get_renormalization_scale)(&mu);
}
inline void get_flavour_scheme_rcl
            (int &nf)
{
  RCLMOD(input,get_flavour_scheme)(&nf);
}
inline void set_alphas_masses_rcl
            (const double mc, const double mb, const double mt,
	           const double wc, const double wb, const double wt)
{
  RCLMOD(input,set_alphas_masses)(&mc,&mb,&mt,&wc,&wb,&wt);
}
inline void set_alphas_masses_rcl
             (const double mc, const double mb, const double mt)
{
  RCLMOD(wrapper,set_alphas_masses_nowidtharg)(&mc,&mb,&mt);
}
inline void use_gfermi_scheme_rcl()
{
  RCLMOD(wrapper,use_gfermi_scheme_noarg)();
}
inline void use_gfermi_scheme_and_set_alpha_rcl
            (const double a)
{
  RCLMOD(wrapper,use_gfermi_scheme_and_set_alpha)(&a);
}
inline void use_gfermi_scheme_and_set_gfermi_rcl
            (const double g)
{
  RCLMOD(wrapper,use_gfermi_scheme_and_set_gfermi)(&g);
}
inline void use_gfermi_scheme_real_and_set_gfermi_rcl
            (const double g)
{
  RCLMOD(wrapper,use_gfermi_scheme_real_and_set_gfermi)(&g);
}
inline void use_gfermi_scheme_complex_and_set_gfermi_rcl
            (const double g)
{
  RCLMOD(wrapper,use_gfermi_scheme_complex_and_set_gfermi)(&g);
}
inline void use_alpha0_scheme_rcl
            (const double a)
{
  RCLMOD(input,use_alpha0_scheme)(&a);
}
inline void use_alpha0_scheme_rcl()
{
  RCLMOD(wrapper,use_alpha0_scheme_noarg)();
}
inline void use_alphaz_scheme_rcl
            (const double a)
{
  RCLMOD(input,use_alphaz_scheme)(&a);
}
inline void use_alphaz_scheme_rcl()
{
  RCLMOD(wrapper,use_alphaz_scheme_noarg)();
}
inline void get_alpha_rcl
            (double &a)
{
  RCLMOD(input,get_alpha)(&a);
}
inline void set_complex_mass_scheme_rcl()
{
  RCLMOD(input,set_complex_mass_scheme)();
}
inline void set_on_shell_scheme_rcl()
{
  RCLMOD(input,set_on_shell_scheme)();
}
inline void set_resonant_particle_rcl
            (const std::string pa)
{
  RCLMOD(input,set_resonant_particle)(pa.c_str(),pa.length());
}
inline void set_quarkline_rcl
            (const int npr, const int q1, const int q2)
{
  RCLMOD(input,set_quarkline)(&npr,&q1,&q2);
}
inline void switchon_resonant_selfenergies_rcl()
{
  RCLMOD(input,switchon_resonant_selfenergies)();
}
inline void switchoff_resonant_selfenergies_rcl()
{
  RCLMOD(input,switchoff_resonant_selfenergies)();
}
inline void set_dynamic_settings_rcl
            (const int n)
{
  RCLMOD(input,set_dynamic_settings)(&n);
}
inline void set_momenta_correction_rcl
            (const bool mc)
{
  int boolint;
  if (mc) {
    boolint = 1;
  } else {
    boolint = 0;
  }
  RCLMOD(input,set_momenta_correction)(&boolint);
}
inline void set_draw_level_branches_rcl
            (const int n)
{
  RCLMOD(input,set_draw_level_branches)(&n);
}
inline void set_print_level_amplitude_rcl
            (const int n)
{
  RCLMOD(input,set_print_level_amplitude)(&n);
}
inline void set_print_level_squared_amplitude_rcl
            (const int n)
{
  RCLMOD(input,set_print_level_squared_amplitude)(&n);
}
inline void set_print_level_correlations_rcl
            (const int n)
{
  RCLMOD(input,set_print_level_correlations)(&n);
}
inline void set_print_level_RAM_rcl
            (const int n)
{
  RCLMOD(input,set_print_level_ram)(&n);
}
inline void scale_coupling3_rcl
            (const std::complex<double> fac,
             const std::string pa1,
             const std::string pa2,
             const std::string pa3)
{
   dcomplex dfac;
   dfac.dr = fac.real();
   dfac.di = fac.imag();
   RCLMOD(input,scale_coupling3)(&dfac,
                                       pa1.c_str(),
                                       pa2.c_str(),
                                       pa3.c_str(),
                                       pa1.length(),
                                       pa2.length(),
                                       pa3.length());
}
inline void scale_coupling4_rcl
            (const std::complex<double> fac,
             const std::string pa1,
             const std::string pa2,
             const std::string pa3,
             const std::string pa4)
{
   dcomplex dfac;
   dfac.dr = fac.real();
   dfac.di = fac.imag();
   RCLMOD(input,scale_coupling4)(&dfac,
                                       pa1.c_str(),
                                       pa2.c_str(),
                                       pa3.c_str(),
                                       pa4.c_str(),
                                       pa1.length(),
                                       pa2.length(),
                                       pa3.length(),
                                       pa4.length());
}
inline void switchoff_coupling3_rcl
            (const std::string pa1,
             const std::string pa2,
             const std::string pa3)
{
  RCLMOD(input,switchoff_coupling3)(pa1.c_str(),
                                          pa2.c_str(),
                                          pa3.c_str(),
                                          pa1.length(),
                                          pa2.length(),
                                          pa3.length());
}
inline void switchoff_coupling4_rcl
            (const std::string pa1,
             const std::string pa2,
             const std::string pa3,
             const std::string pa4)
{
  RCLMOD(input,switchoff_coupling4)(pa1.c_str(),
                                          pa2.c_str(),
                                          pa3.c_str(),
                                          pa4.c_str(),
                                          pa1.length(),
                                          pa2.length(),
                                          pa3.length(),
                                          pa4.length());
}
inline void set_ifail_rcl
            (const int i)
{
   RCLMOD(input,set_ifail)(&i);
}
inline void get_ifail_rcl
            (int &i)
{
   RCLMOD(input,get_ifail)(&i);
}
inline void set_collier_output_dir_rcl
            (const std::string x)
{
  RCLMOD(input,set_collier_output_dir)(x.c_str(),x.length());
}
inline void set_compute_ir_poles_rcl
      (int i)
{
  RCLMOD(input,set_compute_ir_poles)(&i);
}
inline void set_output_file_rcl
            (const std::string x)
{
  RCLMOD(input,set_output_file)(x.c_str(),x.length());
}
/*
  module process_defintion_rcl
*/
inline void define_process_rcl
            (const int npr,
             const std::string process,
             const std::string order)
{
  RCLMOD(process_definition,define_process)(&npr,
                                                  process.c_str(),
                                                  order.c_str(),
                                                  process.length(),
                                                  order.length());
}
inline void set_gs_power_rcl
            (const int npr, const int gsarray[][2], const int gslen)
{
   RCLMOD(wrapper,wrapper_set_gs_power)(&npr,gsarray,&gslen);
}
inline void select_gs_power_BornAmpl_rcl
            (const int npr, const int gspower)
{
  RCLMOD(process_definition,select_gs_power_bornampl)(&npr,&gspower);
}
inline void select_gs_power_LoopAmpl_rcl
            (const int npr, const int gspower)
{
  RCLMOD(process_definition,select_gs_power_loopampl)(&npr,&gspower);
}
inline void unselect_gs_power_BornAmpl_rcl
            (const int npr, const int gspower)
{
  RCLMOD(process_definition,unselect_gs_power_bornampl)(&npr,&gspower);
}
inline void unselect_gs_power_LoopAmpl_rcl
            (const int npr, const int gspower)
{
  RCLMOD(process_definition,unselect_gs_power_loopampl)(&npr,&gspower);
}
inline void select_all_gs_powers_BornAmpl_rcl
            (const int npr)
{
  RCLMOD(process_definition,select_all_gs_powers_bornampl)(&npr);
}
inline void select_all_gs_powers_LoopAmpl_rcl
            (const int npr)
{
  RCLMOD(process_definition,select_all_gs_powers_loopampl)(&npr);
}
inline void unselect_all_gs_powers_BornAmpl_rcl
            (const int npr)
{
  RCLMOD(process_definition,unselect_all_gs_powers_bornampl)(&npr);
}
inline void unselect_all_gs_powers_LoopAmpl_rcl
            (const int npr)
{
  RCLMOD(process_definition,unselect_all_gs_powers_loopampl)(&npr);
}
inline void split_collier_cache_rcl(const int npr, const int n)
{
  RCLMOD(process_definition,split_collier_cache)(&npr,&n);
}
/*
  module process_generation_rcl:
*/
inline void generate_processes_rcl()
{
  RCLMOD(process_generation,generate_processes)();
}
inline void process_exists_rcl(const int npr, bool& exists)
{
  int boolint;
  RCLMOD(process_generation,process_exists)(&npr, &boolint);
  if (boolint) {
    exists = true;
  } else {
    exists = false;
  }
}
/*
  module process_computation_rcl:
*/
inline void set_resonant_squared_momentum_rcl
            (const int npr,const int res, const double ps)
{
  RCLMOD(process_computation,set_resonant_squared_momentum)(&npr,&res,&ps);
}
inline void compute_running_alphas_rcl
            (const double q, const int nf, const int lp)
{
  RCLMOD(process_computation,compute_running_alphas)(&q,&nf,&lp);
}
inline void set_outgoing_momenta_rcl
     (const int npr, double pIn[][4], double p[][4])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  RCLMOD(wrapper,wrapper_set_outgoing_momenta)(&npr, pIn, p, &legs);
}
/*compute_process_rcl is overloaded*/
inline void compute_process_rcl
            (const int npr,
             const double p[][4],
             const std::string order,
             double A2[2],
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_process)(&npr,
                                                p,
                                                &legs,
                                                order.c_str(),
                                                A2,
                                                &boolint,
                                                order.length());
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_process_rcl
            (const int npr,
             const double p[][4],
             const std::string order,
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  double A2[2];
  int boolint;
  RCLMOD(wrapper,wrapper_compute_process)(&npr,
                                                p,
                                                &legs,
                                                order.c_str(),
                                                A2,
                                                &boolint,
                                                order.length());
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_process_rcl
            (const int npr,
             const double p[][4],
             const std::string order,
             double A2[2])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_process)(&npr,
                                                p,
                                                &legs,
                                                order.c_str(),
                                                A2,
                                                &boolint,
                                                order.length());
}
inline void compute_process_rcl
            (const int npr,
             const double p[][4],
             const std::string order)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  double A2[2];
  int boolint;
  RCLMOD(wrapper,wrapper_compute_process)(&npr,
                                                p,
                                                &legs,
                                                order.c_str(),
                                                A2,
                                                &boolint,
                                                order.length());
}
inline void rescale_process_rcl
            (const int npr, const std::string order, double A2[2])
{
  RCLMOD(wrapper,wrapper_rescale_process)(&npr,
                                                order.c_str(),
                                                A2,
                                                order.length());
}
inline void rescale_process_rcl
            (const int npr, const std::string order)
{
  double A2[2];
  RCLMOD(wrapper,wrapper_rescale_process)(&npr,
                                                order.c_str(),
                                                A2,
                                                order.length());
}

inline int get_legs_rcl
            (const int npr)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  return legs;
}

inline int get_cs_rcl
            (const int npr)
{
  int cs;
  RCLMOD(wrapper,get_cs)(&npr, &cs);
  return cs;
}
inline int get_cf_rcl
            (const int npr)
{
  int cf;
  RCLMOD(wrapper,get_cf)(&npr, &cf);
  return cf;
}

inline void get_colour_configurations_rcl
      (int npr, int (*col))
{
  int _legs = get_legs_rcl(npr);
  int _cs = get_cs_rcl(npr);
  RCLMOD(wrapper,wrapper_get_colour_configurations)(&npr,&_legs,&_cs,col);
}

inline void get_helicity_configurations_rcl
      (int npr, int (*col))
{
  int _legs = get_legs_rcl(npr);
  int _cf = get_cf_rcl(npr);
  RCLMOD(wrapper,wrapper_get_helicity_configurations)(&npr,&_legs,&_cf,col);
}

inline void get_amplitude_rcl
            (const int npr,
             const int pow,
             const std::string order,
             const int colour[], const int hel[],
             std::complex<double>& A)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  dcomplex dA;
  dA.di = A.imag();
  dA.dr = A.real();
  RCLMOD(wrapper,wrapper_get_amplitude)(&npr,
                                              &pow,
                                              order.c_str(),
                                              colour,
                                              hel,
                                              &legs,
                                              &dA,
                                              order.length());
  A = std::complex<double> (dA.dr, dA.di);
}
inline void get_squared_amplitude_rcl
            (const int npr,
             const int pow,
             const std::string order,
             double &A2)
{
  RCLMOD(process_computation,get_squared_amplitude)(&npr,
                                                          &pow,
                                                          order.c_str(),
                                                          &A2,
                                                          order.length());
}
inline void get_polarized_squared_amplitude_rcl
            (const int npr,
             const int pow,
             const std::string order,
             const int hel[],
             double &A2h)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  RCLMOD(wrapper,wrapper_get_polarized_squared_amplitude)(&npr,
                                                                &pow,
                                                                order.c_str(),
                                                                hel,
                                                                &legs,
                                                                &A2h,
                                                                order.length());
}
inline void compute_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             double &A2cc,
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_colour_correlation)(&npr,
                                                           p,
                                                           &legs,
                                                           &i1,
                                                           &i2,
                                                           &A2cc,
                                                           &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2cc;
  RCLMOD(wrapper,wrapper_compute_colour_correlation)(&npr,
                                                           p,
                                                           &legs,
                                                           &i1,
                                                           &i2,
                                                           &A2cc,
                                                           &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             double &A2cc)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_colour_correlation)(&npr,
                                                           p,
                                                           &legs,
                                                           &i1,
                                                           &i2,
                                                           &A2cc,
                                                           &boolint);
}
inline void compute_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2cc;
  RCLMOD(wrapper,wrapper_compute_colour_correlation)(&npr,
                                                           p,
                                                           &legs,
                                                           &i1,
                                                           &i2,
                                                           &A2cc,
                                                           &boolint);
}
inline void compute_colour_correlation_int_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             double &A2ccint,
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_colour_correlation_int)(&npr,
                                                               p,
                                                               &legs,
                                                               &i1,
                                                               &i2,
                                                               &A2ccint,
                                                               &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_colour_correlation_int_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2ccint;
  RCLMOD(wrapper,wrapper_compute_colour_correlation_int)(&npr,
                                                               p,
                                                               &legs,
                                                               &i1,
                                                               &i2,
                                                               &A2ccint,
                                                               &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_colour_correlation_int_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             double &A2ccint)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_colour_correlation_int)(&npr,
                                                               p,
                                                               &legs,
                                                               &i1,
                                                               &i2,
                                                               &A2ccint,
                                                               &boolint);
}
inline void compute_colour_correlation_int_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2ccint;
  RCLMOD(wrapper,wrapper_compute_colour_correlation_int)(&npr,
                                                               p,
                                                               &legs,
                                                               &i1,
                                                               &i2,
                                                               &A2ccint,
                                                               &boolint);
}
inline void compute_all_colour_correlations_rcl
            (const int npr,
             const double p[][4],
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_all_colour_correlations)(&npr,
                                                                p,
                                                                &legs,
                                                                &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_all_colour_correlations_rcl
            (const int npr,
             const double p[][4])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_all_colour_correlations)(&npr,
                                                                p,
                                                                &legs,
                                                                &boolint);
}
inline void compute_all_colour_correlations_int_rcl
            (const int npr,
             const double p[][4],
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_all_colour_correlations_int)(&npr,
                                                                    p,
                                                                    &legs,
                                                                    &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_all_colour_correlations_int_rcl
            (const int npr,
             const double p[][4])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_all_colour_correlations_int)(&npr,
                                                                    p,
                                                                    &legs,
                                                                    &boolint);
}
inline void rescale_colour_correlation_rcl
            (const int npr,
             const int i1, const int i2,
             double &A2cc)
{
  RCLMOD(wrapper,wrapper_rescale_colour_correlation)(&npr,
                                                           &i1,
                                                           &i2,
                                                           &A2cc);
}
inline void rescale_colour_correlation_rcl
            (const int npr,
             const int i1, const int i2)
{
  double A2cc;
  RCLMOD(wrapper,wrapper_rescale_colour_correlation)(&npr,
                                                           &i1,
                                                           &i2,
                                                           &A2cc);
}
inline void rescale_colour_correlation_int_rcl
            (const int npr,
             const int i1, const int i2,
             double &A2ccint)
{
  RCLMOD(wrapper,wrapper_rescale_colour_correlation_int)(&npr,
                                                               &i1,
                                                               &i2,
                                                               &A2ccint);
}
inline void rescale_colour_correlation_int_rcl
            (const int npr,
             const int i1, const int i2)
{
  double A2ccint;
  RCLMOD(wrapper,wrapper_rescale_colour_correlation_int)(&npr,
                                                               &i1,
                                                               &i2,
                                                               &A2ccint);
}
inline void rescale_all_colour_correlations_rcl
            (const int npr)
{
  RCLMOD(process_computation,rescale_all_colour_correlations)(&npr);
}
inline void get_colour_correlation_rcl
            (const int npr,
             const int pow,
             const int i1, const int i2,
             double& A2cc)
{
  RCLMOD(process_computation,get_colour_correlation)(&npr,
                                                           &pow,
                                                           &i1,
                                                           &i2,
                                                           &A2cc);
}
inline void get_colour_correlation_int_rcl
            (const int npr,
             const int pow,
             const int i1, const int i2,
             double& A2ccint)
{
  RCLMOD(process_computation,get_colour_correlation_int)(&npr,
                                                               &pow,
                                                               &i1,
                                                               &i2,
                                                               &A2ccint);
}
inline void compute_spin_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             const std::complex<double> v[4],
             double &A2scc,
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  int boolint;
  RCLMOD(wrapper,wrapper_compute_spin_colour_correlation)(&npr,
                                                                p,
                                                                &legs,
                                                                &i1,
                                                                &i2,
                                                                dfac,
                                                                &A2scc,
                                                                &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_spin_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             const std::complex<double> v[4],
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  double A2scc;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  int boolint;
  RCLMOD(wrapper,wrapper_compute_spin_colour_correlation)(&npr,
                                                                p,
                                                                &legs,
                                                                &i1,
                                                                &i2,
                                                                dfac,
                                                                &A2scc,
                                                                &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_spin_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             const std::complex<double> v[4],
             double &A2scc)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_compute_spin_colour_correlation)(&npr,
                                                                p,
                                                                &legs,
                                                                &i1,
                                                                &i2,
                                                                dfac,
                                                                &A2scc,
                                                                &boolint);
}
inline void compute_spin_colour_correlation_rcl
            (const int npr,
             const double p[][4],
             const int i1, const int i2,
             const std::complex<double> v[4])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2scc;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_compute_spin_colour_correlation)(&npr,
                                                                p,
                                                                &legs,
                                                                &i1,
                                                                &i2,
                                                                dfac,
                                                                &A2scc,
                                                                &boolint);
}
inline void rescale_spin_colour_correlation_rcl
            (const int npr,
             const int i1, const int i2,
             const std::complex<double> v[4],
             double &A2scc)
{
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_rescale_spin_colour_correlation)(&npr,
                                                                &i1,
                                                                &i2,
                                                                dfac,
                                                                &A2scc);
}
inline void rescale_spin_colour_correlation_rcl
            (const int npr,
             const int i1, const int i2,
             const std::complex<double> v[4])
{
  double A2scc;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_rescale_spin_colour_correlation)(&npr,
                                                                &i1,
                                                                &i2,
                                                                dfac,
                                                                &A2scc);
}
inline void get_spin_colour_correlation_rcl
            (const int npr,
             const int pow,
             const int i1, const int i2,
             double &A2scc)
{
  RCLMOD(process_computation,get_spin_colour_correlation)(&npr,
                                                                &pow,
                                                                &i1,
                                                                &i2,
                                                                &A2scc);
}
inline void compute_spin_correlation_rcl
            (const int npr,
             const double p[][4],
             const int j,
             const std::complex<double> v[4],
             double &A2sc,
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  int boolint;
  RCLMOD(wrapper,wrapper_compute_spin_correlation)(&npr,
                                                         p,
                                                         &legs,
                                                         &j,
                                                         dfac,
                                                         &A2sc,
                                                         &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_spin_correlation_rcl
            (const int npr,
             const double p[][4],
             const int j,
             const std::complex<double> v[4],
             bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  double A2sc;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  int boolint;
  RCLMOD(wrapper,wrapper_compute_spin_correlation)(&npr,
                                                         p,
                                                         &legs,
                                                         &j,
                                                         dfac,
                                                         &A2sc,
                                                         &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_spin_correlation_rcl
            (const int npr,
             const double p[][4],
             const int j,
             const std::complex<double> v[4],
             double &A2sc)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_compute_spin_correlation)(&npr,
                                                         p,
                                                         &legs,
                                                         &j,
                                                         dfac,
                                                         &A2sc,
                                                         &boolint);
}
inline void compute_spin_correlation_rcl
            (const int npr,
             const double p[][4],
             const int j,
             const std::complex<double> v[4])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2sc;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_compute_spin_correlation)(&npr,
                                                         p,
                                                         &legs,
                                                         &j,
                                                         dfac,
                                                         &A2sc,
                                                         &boolint);
}
inline void rescale_spin_correlation_rcl
            (const int npr,
             const int j,
             const std::complex<double> v[4],
             double &A2sc)
{
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_rescale_spin_correlation)(&npr,
                                                         &j,
                                                         dfac,
                                                         &A2sc);
}
inline void rescale_spin_correlation_rcl
            (const int npr,
             const int j,
             const std::complex<double> v[4])
{
  double A2sc;
  dcomplex dfac[4];
  for (int i = 0; i < 4; ++i) {
    dfac[i].dr = v[i].real();
    dfac[i].di = v[i].imag();
  }
  RCLMOD(wrapper,wrapper_rescale_spin_correlation)(&npr,
                                                         &j,
                                                         dfac,
                                                         &A2sc);
}
inline void get_spin_correlation_rcl
            (const int npr, const int pow, double &A2sc)
{
  RCLMOD(process_computation,get_spin_correlation)(&npr,
                                                         &pow,
                                                         &A2sc);
}
inline void compute_spin_correlation_matrix_rcl
      (const int npr,
       const double p[][4],
       const int j,
       double A2scm[4][4],
       bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_spin_correlation_matrix)(&npr,
                                                          p,
                                                          &legs,
                                                          &j,
                                                          A2scm,
                                                          &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_spin_correlation_matrix_rcl
      (const int npr,
       const double p[][4],
       const int j,
       double A2scm[4][4])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  RCLMOD(wrapper,wrapper_compute_spin_correlation_matrix)(&npr,
                                                          p,
                                                          &legs,
                                                          &j,
                                                          A2scm,
                                                          &boolint);
}
inline void compute_spin_correlation_matrix_rcl
      (const int npr,
       const double p[][4],
       const int j,
       bool& momenta_check)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2scm[4][4];
  RCLMOD(wrapper,wrapper_compute_spin_correlation_matrix)(&npr,
                                                          p,
                                                          &legs,
                                                          &j,
                                                          A2scm,
                                                          &boolint);
  if (boolint) {
    momenta_check = true;
  } else {
    momenta_check = false;
  }
}
inline void compute_spin_correlation_matrix_rcl
      (const int npr,
       const double p[][4],
       const int j)
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  int boolint;
  double A2scm[4][4];
  RCLMOD(wrapper,wrapper_compute_spin_correlation_matrix)(&npr,
                                                          p,
                                                          &legs,
                                                          &j,
                                                          A2scm,
                                                          &boolint);
}
inline void rescale_spin_correlation_matrix_rcl
      (const int npr,
       const int j,
       double A2scm[4][4])
{
  RCLMOD(wrapper,wrapper_rescale_spin_correlation_matrix)(&npr,
                                                          &j,
                                                          A2scm);
}
inline void rescale_spin_correlation_matrix_rcl
      (const int npr,
       const int j)
{
  double A2scm[4][4];
  RCLMOD(wrapper,wrapper_rescale_spin_correlation_matrix)(&npr,
                                                          &j,
                                                          A2scm);
}
inline void get_spin_correlation_matrix_rcl
      (const int npr,
       int pow,
       double A2scm[4][4])
{
  RCLMOD(process_computation,get_spin_correlation_matrix)(&npr,
                                                          &pow,
                                                          A2scm);
}
inline void get_momenta_rcl
            (const int npr, double p[][4])
{
  int legs;
  RCLMOD(wrapper,get_legs)(&npr, &legs);
  RCLMOD(wrapper,wrapper_get_momenta)(&npr,p,&legs);
}
inline void get_recola_version_rcl
      (std::string& recolaversion)
{
  char rv_cstr[10];
  int len;
  RCLMOD(wrapper,wrapper_get_recola_version)(rv_cstr,&len,10);
  recolaversion = std::string(rv_cstr,len);
}
inline void set_tis_required_accuracy_rcl
            (const double acc)
{
  RCLMOD(process_computation,set_tis_required_accuracy)(&acc);
}
inline void get_tis_required_accuracy_rcl
            (double &acc)
{
  RCLMOD(process_computation,get_tis_required_accuracy)(&acc);
}
inline void set_tis_critical_accuracy_rcl
            (const double acc)
{
  RCLMOD(process_computation,set_tis_critical_accuracy)(&acc);
}
inline void get_tis_critical_accuracy_rcl
            (double &acc)
{
  RCLMOD(process_computation,get_tis_critical_accuracy)(&acc);
}
inline void get_tis_accuracy_flag_rcl
            (int &flag)
{
  RCLMOD(process_computation,get_tis_accuracy_flag)(&flag);
}
/*
  module reset_rcl
*/

inline void reset_recola_rcl()
{
  RCLMOD(reset,reset_recola)();
}

}

#endif // RECOLA_H_
