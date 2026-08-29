import SubstrateFramework

/-! The infrastructure theorem must stay axiom-free; each ingested
    file's main (lead) theorem is audited alongside with its full
    axiom footprint (the standard Mathlib trio is expected). -/

#print axioms SubstrateFramework.compose_implications
#print axioms SubstrateFramework.OpticalGothic.commonMetricComposition

#print axioms ActionQuantum.sqrt_one_sub_sq_pos  -- ActionQuantum.lean
#print axioms SGAmplitudeCondition.amplitudeCondition_feasible  -- AmplitudeCondition.lean
#print axioms SGBasin.basin_single_point_exists  -- Basin.lean
#print axioms SGBasin.basin_positive_measure  -- Basin.lean (load-bearing promoted interval theorem)
#print axioms SGBridge.exchange_phase_free_fermion  -- Bridge.lean
#print axioms SGChargeDiscrimination.channels_cp_conjugate  -- ChargeDiscrimination.lean
#print axioms ComparsiVirial.periodic_virial_forces_zero  -- ComparsiVirial.lean
#print axioms SGDetectorGeometry.kinkComptonTime_eq_inv_mass  -- DetectorGeometry.lean
#print axioms SGDoubleWell.barrier_height  -- DoubleWell.lean
#print axioms SGDrivenNeumannBC.drivenNeumannBC_periodic  -- DrivenNeumannBC.lean
#print axioms SGEnergy.binding_gap_nonneg  -- Energy.lean
#print axioms SGEnergyFlux.energyFlux_sign_absorbing  -- EnergyFlux.lean
#print axioms SGExchange.kink_exchange_neg  -- Exchange.lean
#print axioms SGFormalization.binding_gap_nonneg  -- Formalization.lean
#print axioms SGGates.gate1_energy_budget  -- Gates.lean
#print axioms SGHBTApproach.c_cre_mul_c_ann  -- HBT_Approach.lean
#print axioms SGLandauZener.landauZener_prob_range  -- LandauZener.lean
#print axioms Phase10CF.center_card_three  -- Phase10CF_Z3CenterConfinement.lean
#print axioms Phase11Flavor.cpPhases_two  -- Phase11Flavor_CPPhaseCount.lean
#print axioms Phase12GW.gravity_lowest_multipole_two  -- Phase12GW_LowestMultipole.lean
#print axioms Phase13FS.quadratic_source_radiates_at_2omega  -- Phase13FS_HarmonicDoubling.lean
#print axioms Phase14P3D.spherical_breather_radiates_nothing  -- Phase14P3D_SphericalNull.lean
#print axioms Phase15NC.topo_conserved  -- Phase15NC_NonlinearChiralSplit.lean
#print axioms Phase16QB.axisym_one_polarization  -- Phase16QB_TwoPolarizations.lean
#print axioms Phase17WZW.wzw_level_eq_baryon_winding  -- Phase17WZW_LevelWinding.lean
#print axioms Phase18Chiral.pion_triplet  -- Phase18Chiral_GoldstoneCount.lean
#print axioms Phase19D3S.analytic_even_symbol_power_two_is_s1  -- Phase19_D3S_LocalityImpliesS1.lean
#print axioms Phase19OM1.chiF_generator  -- Phase19_OM1_Z2Character.lean
#print axioms Phase1ThermalGate.thermal_gate_eq_sech_sq  -- Phase1Bridge_ThermalGate.lean
#print axioms Phase20ME.halfQuantum_add_self  -- Phase20ME_HalfQuantumOrder2.lean
#print axioms Phase20MH.overlap_pos  -- Phase20MH_OverlapHierarchy.lean
#print axioms Phase21AS.odV1_underdetermined  -- Phase21AS_OverDetermination.lean
#print axioms Phase22AS.odV2_overdetermined  -- Phase22AS_BetaPinnedScale.lean
#print axioms Phase23EW.tr_T3sq_eq  -- Phase23EW_Sin2ThetaW.lean
#print axioms Phase25Seeding.seeding_dc_null  -- Phase25SeedingKernel.lean
#print axioms Phase26Lifetime.lifetime_strictAntiOn  -- Phase26LifetimeKernel.lean
#print axioms Phase27Medium.breather_exists_example  -- Phase27MediumKernel.lean
#print axioms Phase28.preRungC_all_static  -- Phase28BarrierKernel.lean
#print axioms Phase29.bps_zero_binding  -- Phase29YieldKernel.lean
#print axioms Phase30PN.subdivision_bound  -- Phase30PNKernel.lean
#print axioms Phase31CM.collective_exceeds_two_body  -- Phase31CMKernel.lean
#print axioms Phase32GB.enh_grows  -- Phase32GBKernel.lean
#print axioms Phase33.derivedScale_flow_invariant  -- Phase34KIKernel.lean
#print axioms Phase37WN.weight_pos  -- Phase37WNKernel.lean
#print axioms Phase38MD.overparam_cancel  -- Phase38MDKernel.lean
#print axioms Phase3EM.u1_divergence_identity  -- Phase3EM_U1Current.lean
#print axioms Phase40.quad_traceless  -- Phase40TX_RotatingTorus.lean
#print axioms Phase41.induced_normalization_closes  -- Phase41GaugeKinetic.lean
#print axioms Phase46Lepton.parity_eq_ite  -- Phase46EL_LeptonConstraint.lean
#print axioms Phase47BM.gap_eq  -- Phase47BM_RigidTopSpectrum.lean
#print axioms Phase48CE.gap_independent_of_vacuum_energy  -- Phase48CE_OperatorAndGap.lean
#print axioms Phase4Strong.skyrmion_odd_is_fermion  -- Phase4Strong_FRSpinStat.lean
#print axioms Phase5Gravity.graviton_D4_eq_two  -- Phase5Gravity_GravitonTT.lean
#print axioms Phase6Weak.parity_involution  -- Phase6Weak_MaxParityViolation.lean
#print axioms Phase7EW.Qop_eq_T3_plus_Yhalf  -- Phase7EW_MasslessPhoton.lean
#print axioms Phase8QCD.adjoint_dim_eq  -- Phase8QCD_SU3GaugeFacts.lean
#print axioms Phase9SM.anomaly_U1_cubed  -- Phase9SM_AnomalyCancellation.lean
#print axioms SGPhaseCondition.optimal_phase_satisfies_condition  -- PhaseCondition.lean
#print axioms SGProductionAmplitude.sMatrixTimeDelay_free_fermion  -- ProductionAmplitude.lean
#print axioms SGSigma2Prod.sigma2_total_bounded  -- Sigma2Prod.lean
