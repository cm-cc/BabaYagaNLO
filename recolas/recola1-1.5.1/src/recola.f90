!#####################################################################
!!
!!  File  recola.f90
!!  is part of RECOLA (REcursive Computation of One Loop Amplitudes)
!!
!!  Copyright (C) 2015-2025   Stefano Actis, Ansgar Denner,
!!                            Lars Hofer, Jean-Nicolas Lang,
!!                            Andreas Scharf, Sandro Uccirati
!!
!!  RECOLA is licenced under the GNU GPL version 3,
!!         see COPYING for details.
!!
!#####################################################################

  module recola

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use globals_rcl, only: &
      set_pure_QED_rcl, &
      unset_pure_QED_rcl, &
      get_recola_version_rcl, &
      print_warning_summary_rcl

  use input_rcl, only : &
      set_pole_mass_z_rcl,set_onshell_mass_z_rcl, &
      set_pole_mass_w_rcl,set_onshell_mass_w_rcl, &
      set_pole_mass_h_rcl, &
      set_pole_mass_electron_rcl, &
      set_pole_mass_muon_rcl, &
      set_pole_mass_tau_rcl, &
      set_pole_mass_up_rcl, &
      set_pole_mass_down_rcl, &
      set_pole_mass_charm_rcl, &
      set_pole_mass_strange_rcl, &
      set_pole_mass_top_rcl, &
      set_pole_mass_bottom_rcl, &
      set_light_fermions_rcl, &
      set_light_electron_rcl, &
      set_light_muon_rcl, &
      set_light_tau_rcl, &
      set_light_up_rcl, &
      set_light_down_rcl, &
      set_light_charm_rcl, &
      set_light_strange_rcl, &
      set_light_top_rcl, &
      set_light_bottom_rcl, &
      unset_light_electron_rcl, &
      unset_light_muon_rcl, &
      unset_light_tau_rcl, &
      unset_light_up_rcl, &
      unset_light_down_rcl, &
      unset_light_charm_rcl, &
      unset_light_strange_rcl, &
      unset_light_top_rcl, &
      unset_light_bottom_rcl, &
      use_dim_reg_soft_rcl, &
      use_mass_reg_soft_rcl, &
      set_mass_reg_soft_rcl, &
      set_delta_uv_rcl, &
      set_mu_uv_rcl, &
      set_delta_ir_rcl, &
      set_mu_ir_rcl, &
      set_complex_mass_scheme_rcl, &
      set_on_shell_scheme_rcl, &
      set_alphas_rcl, &
      get_alphas_rcl, &
      get_renormalization_scale_rcl, &
      get_flavour_scheme_rcl, &
      set_alphas_masses_rcl, &
      use_gfermi_scheme_rcl, &
      use_alpha0_scheme_rcl, &
      use_alphaz_scheme_rcl, &
      use_alphaMSbar_scheme_rcl, &
      set_gfermi_rcl, &
      set_alpha0_rcl, &
      set_alphaz_rcl, &
      set_alphaMSbar_rcl, &
      get_alpha_rcl, &
      set_resonant_particle_rcl, &
      switchon_resonant_selfenergies_rcl, &
      switchoff_resonant_selfenergies_rcl, &
      switchon_transverse_resonant_propagator_rcl, &
      switchoff_transverse_resonant_propagator_rcl, &
      set_dynamic_settings_rcl, &
      set_momenta_correction_rcl, &
      set_draw_level_branches_rcl, &
      set_print_level_amplitude_rcl, &
      set_print_level_squared_amplitude_rcl, &
      set_print_level_correlations_rcl, &
      set_print_level_RAM_rcl, &
      scale_coupling3_rcl, &
      scale_coupling4_rcl, &
      switchoff_coupling3_rcl, &
      switchoff_coupling4_rcl, &
      set_ifail_rcl, get_ifail_rcl, &
      set_collier_mode_rcl, &
      set_collier_output_dir_rcl, &
      set_compute_ir_poles_rcl, &
      set_output_file_rcl

  use process_definition_rcl, only : &
      define_process_rcl, &
      set_offshell_photons_rcl, &
      use_gfermi_mixed_scheme_rcl, &
      use_alpha0_mixed_scheme_rcl, &
      use_alphaz_mixed_scheme_rcl, &
      use_alphaMSbar_mixed_scheme_rcl, &
      get_reference_alpha_rcl, &
      get_NLOalpha_rcl, &
      get_alphaGF_rcl, &
      get_alpha0_rcl, &
      get_alphaz_rcl, &
      get_alphaMSbar_rcl, &
      get_QrenMSbar_rcl, &
      set_gs_power_rcl, &
      select_gs_power_BornAmpl_rcl, &
      select_gs_power_LoopAmpl_rcl, &
      unselect_gs_power_BornAmpl_rcl, &
      unselect_gs_power_LoopAmpl_rcl, &
      select_all_gs_powers_BornAmpl_rcl, &
      select_all_gs_powers_LoopAmpl_rcl, &
      unselect_all_gs_powers_BornAmpl_rcl, &
      unselect_all_gs_powers_LoopAmpl_rcl, &
      set_internal_projection_rcl, &
      set_quarkline_rcl, &
      split_collier_cache_rcl

  use process_generation_rcl, only : &
      generate_processes_rcl,        &
      get_deltar_rcl, &
      get_dalZ_rcl, &
      get_dalMS_rcl, &
      process_exists_rcl

  use process_computation_rcl, only : &
      set_resonant_squared_momentum_rcl, &
      compute_running_alphas_rcl, &
      compute_process_rcl, &
      rescale_process_rcl, &
      get_amplitude_rcl, &
      get_squared_amplitude_rcl, &
      get_polarized_squared_amplitude_rcl, &
      compute_colour_correlation_rcl, &
      rescale_colour_correlation_rcl, &
      compute_all_colour_correlations_rcl, &
      rescale_all_colour_correlations_rcl, &
      get_colour_correlation_rcl, &
      compute_spin_colour_correlation_rcl, &
      rescale_spin_colour_correlation_rcl, &
      get_spin_colour_correlation_rcl, &
      compute_spin_correlation_rcl, &
      rescale_spin_correlation_rcl, &
      get_spin_correlation_rcl, &
      compute_spin_correlation_matrix_rcl, &
      rescale_spin_correlation_matrix_rcl, &
      get_spin_correlation_matrix_rcl, &
      compute_int_colour_correlation_rcl, &
      rescale_int_colour_correlation_rcl, &
      compute_all_int_colour_correlations_rcl, &
      rescale_all_int_colour_correlations_rcl, &
      get_int_colour_correlation_rcl, &
      compute_nlo_colour_correlation_rcl, &
      rescale_nlo_colour_correlation_rcl, &
      get_nlo_colour_correlation_rcl, &
      get_momenta_rcl, &
      get_colour_configurations_rcl, &
      get_helicity_configurations_rcl, &
      set_TIs_required_accuracy_rcl, &
      get_TIs_required_accuracy_rcl, &
      set_TIs_critical_accuracy_rcl, &
      get_TIs_critical_accuracy_rcl, &
      get_TIs_accuracy_flag_rcl

  use outgoing_momenta_rcl, only: &
      set_outgoing_momenta_rcl

  use reset_rcl, only : &
      reset_recola_rcl

  use collier

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module recola

!#########################################################################
