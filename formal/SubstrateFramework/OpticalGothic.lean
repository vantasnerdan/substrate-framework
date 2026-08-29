import Mathlib

/-!
Finite algebraic and logical glue for the P248 optical-gothic audit.

The block-matrix bijection and its ten-component Jacobian are proved primarily
by the campaign's exact SymPy verifier.  This file corroborates the scalar
compatibility, sign, shear, and theorem-composition steps without encoding an
ether ontology or a microscopic constitutive action.
-/

namespace SubstrateFramework.OpticalGothic

def materialResidual (n nt adv div : ℝ) : ℝ := nt + adv + n * div

def harmonicResidual (s n nt adv div : ℝ) : ℝ :=
  2 * n * nt - s * (2 * n * adv + n ^ 2 * div)

/-- With the material-flow shift sign, gothic continuity differs from twice
material continuity by the compression term. -/
theorem materialSignCompatibility (n nt adv div : ℝ) :
    harmonicResidual (-1) n nt adv div =
      2 * n * materialResidual n nt adv div - n ^ 2 * div := by
  simp [harmonicResidual, materialResidual]
  ring

/-- For positive determinant mean, simultaneous material and harmonic
continuity in the material-flow convention forces incompressibility. -/
theorem jointContinuityForcesIncompressible
    (n nt adv div : ℝ) (hn : 0 < n)
    (hm : materialResidual n nt adv div = 0)
    (hh : harmonicResidual (-1) n nt adv div = 0) : div = 0 := by
  rw [materialSignCompatibility, hm] at hh
  have hn2 : 0 < n ^ 2 := sq_pos_of_pos hn
  nlinarith

/-- Spatial inhomogeneity is not excluded: setting the local temporal,
advective, and compressive derivatives to zero satisfies both residuals for
every positive local value of the optical determinant mean. -/
theorem staticInhomogeneousWitness (n : ℝ) :
    materialResidual n 0 0 0 = 0 ∧ harmonicResidual (-1) n 0 0 0 = 0 := by
  simp [materialResidual, harmonicResidual]

/-- The plus shift sign printed in the source gives a different joint
compatibility condition when material continuity is imposed. -/
theorem printedSignCompatibility (n adv div : ℝ) :
    harmonicResidual 1 n (-adv - n * div) adv div =
      -n * (4 * adv + 3 * n * div) := by
  simp [harmonicResidual]
  ring

/-- The proposed additive determinant-mean strain has a nonzero trace on the
positive diagonal witness (4,2,1), whose determinant mean is 2. -/
theorem additiveShearCounterexample :
    (4 - 2 : ℝ) + (2 - 2) + (1 - 2) = 1 := by norm_num

/-- Subtracting the mean logarithm makes the three spectral shear coordinates
exactly trace-free. -/
theorem logarithmicShearTrace (l₁ l₂ l₃ : ℝ) :
    let mean := (l₁ + l₂ + l₃) / 3
    (l₁ - mean) + (l₂ - mean) + (l₃ - mean) = 0 := by
  dsimp
  ring

/-- A strictly positive magnitude cannot simultaneously be the positive and
negative Newtonian field-energy sign. -/
theorem newtonianEnergySignConflict (energy : ℝ) (henergy : 0 < energy) :
    energy ≠ -energy := by
  nlinarith

/-- The complete point-map Jacobian factor is nonzero for positive lapse and
positive spatial determinant. -/
theorem completeJacobianNonzero (lapse spatialDet : ℝ)
    (hlapse : 0 < lapse) (hdet : 0 < spatialDet) :
    2 * lapse * spatialDet ≠ 0 := by
  positivity

/-- Reusing the on-shell self-source as an arbitrary-field definition cancels
the field operator and leaves only the matter source. -/
theorem onShellClosureCancellation (boxField coupling matter : ℝ)
    (hcoupling : coupling ≠ 0) :
    boxField + coupling * (matter - boxField / coupling) = coupling * matter := by
  field_simp [hcoupling]
  ring

/-- Logical composition at the fixed synthesis boundary: an invertible
variational pullback transfers the metric Euler equation, while the accepted
massive square-root/Hamiltonian and massless null/affine supplied-metric
sectors retain their exact scopes. -/
theorem commonMetricComposition
    {MetricEuler OpticalEuler MassiveSector MasslessNullAffineSector : Prop}
    (pullback : OpticalEuler ↔ MetricEuler)
    (metricFieldEquation : MetricEuler)
    (massivePremise : MassiveSector)
    (masslessPremise : MasslessNullAffineSector) :
    OpticalEuler ∧ MassiveSector ∧ MasslessNullAffineSector := by
  exact ⟨pullback.mpr metricFieldEquation,
    massivePremise,
    masslessPremise⟩

end SubstrateFramework.OpticalGothic
