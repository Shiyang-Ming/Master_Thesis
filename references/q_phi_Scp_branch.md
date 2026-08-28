# Purple S_CP(pi0 K_S) constraint for q-phi plot

## Goal

We want to add the purple branch/band in the q-phi plane.

The purple constraint comes from the mixing-induced CP asymmetry

S_CP(B_d -> pi0 K_S)

and is converted into a constraint on the electroweak penguin parameters:

q and phi.

The logical chain is:

S_CP(pi0 K_S)
  -> phi00
  -> q(phi)
  -> purple branch/band.

This is a phase constraint, not a branching-ratio constraint.

---

## Project variable naming

Use the project variable names consistently.

### Experimental CP observables

Paper notation:

S_CP(pi0 K_S)

Project variable:

Scppi0Ks0
Scppi0KsErr

Paper notation:

A_CP(pi0 K_S)

Project variable:

Acppi0Ks0
Acppi0KsErr

For the old reference figure, the label says:

S_CP(pi0 K_S) = 0.58

However, the current project input may use a newer value, for example:

Scppi0Ks0 = 0.64
Scppi0KsErr = 0.13

When reproducing the historical reference figure, allow an explicit option such as:

Scppi0KsCentralForPlot = 0.58

but do not silently overwrite the global input Scppi0Ks0.

---

### CKM and mixing phases

Paper notation:

gamma

Project variable:

gamma0
gammaErr

Paper notation:

phi_d

Project variable:

phid0
phidErr

All angles are stored internally in radians using Mathematica's Degree multiplier.

---

### B -> pi K hadronic parameters

Paper notation:

r

Project variable:

r0
rErr

Paper notation:

delta

Project variable:

delta0
deltaErr

Paper notation:

r_c

Project variable:

rc0
rcErr

Paper notation:

delta_c

Project variable:

deltac0
deltacErr

---

### Electroweak penguin parameters

Paper notation:

q

Project variable:

q0
qErr

Paper notation:

phi

Project variable:

phi0
phiErr

In the q-phi scan, phi is the horizontal scan variable.

Use phiScan or phiNP as the local scan variable to avoid confusing it with the global input phi0.

The x-axis is phi in degrees.

Internally use radians.

---

## Relation between S_CP(pi0 K_S) and phi00

Use the convention already used in the project if available.

The usual structure is

S_CP(pi0 K_S) =
sqrt(1 - A_CP(pi0 K_S)^2) Sin[phi_d - phi00]

or an equivalent convention depending on the sign conventions of the paper/project.

Important:

Do not assume a single phi00 solution at first.

For a given S_CP, A_CP and phi_d, generate all allowed phi00 branches.

Only after solving for q(phi), apply physical cuts and identify which branch contributes to the visible purple region.

Physical cuts:

- q real
- q > 0
- q inside plot range
- phi inside plot range
- numerical solution valid

---

## Paper formula for the purple constraint

The paper writes

q = (-B_c + Sqrt[B_c^2 - 4 A_c D_c])/(2 A_c)

where

A_c = r_c^2 (-Tan[phi00] Cos[2 phi] - Sin[2 phi])

B_c =
  2 r_c Cos[delta_c] (Tan[phi00] Cos[phi] + Sin[phi])
  - (4/3) chat_plus A_c
  - (2 r_c^2 - 2 r_c r Cos[delta_c - delta])
    (-Tan[phi00] Cos[gamma + phi] - Sin[gamma + phi])

D_c =
  -Tan[phi00]
  - (2 r_c Cos[delta_c] - 2 r Cos[delta])
    (Tan[phi00] Cos[gamma] + Sin[gamma])
  + (r_c^2 + r^2 - 2 r_c r Cos[delta_c - delta])
    (-Tan[phi00] Cos[2 gamma] - Sin[2 gamma])
  + (4/3) atildeC q r_c
    (-Tan[phi00] Cos[phi] - Sin[phi])
  + (4/9) q^2 (atildeS^2 + atildeC^2) A_c
  + (4/3) (-Tan[phi00] Cos[gamma + phi] - Sin[gamma + phi])
    (r_c^2 chat_plus - r_c r (atildeC Cos[delta] + atildeS Sin[delta]))

with

chat_plus = atildeC q Cos[delta_c] + atildeS q Sin[delta_c].

---

## Important notation warning

The paper's A_c is a quadratic-equation coefficient.

Do NOT name it Ac in the code, because the project already uses names like Acppi0Ks0 for CP asymmetries.

Use clear internal names such as:

quadAc
quadBc
quadDc

or

coeffAc
coeffBc
coeffDc

Similarly, avoid naming tilde a_C as Ac.

Use names such as:

atildeC
atildeS

or, for central values,

atildeC0
atildeS0

---

## Definition of atildeC and atildeS used in this project

In this implementation, atildeC and atildeS are not treated as independent input parameters.

They should be extracted from the approximate first-order expressions for R and A_CP(B_d -> pi- K+).

The relevant formulas are:

R =
  1
  - 2 r Cos[delta] Cos[gamma]
  + 2 r_c atildeC q Cos[phi]
  - 2 rho_c Cos[theta_c] Cos[gamma]

A_CP(pi- K+) =
  (4/3) r_c atildeS q Sin[phi]
  - 2 r Sin[delta] Sin[gamma]

where higher-order terms are neglected, following the paper formula.

Therefore,

atildeC =
  (R - 1 + 2 r Cos[delta] Cos[gamma]
     + 2 rho_c Cos[theta_c] Cos[gamma])
  /(2 r_c q Cos[phi])

and

atildeS =
  (A_CP(pi- K+) + 2 r Sin[delta] Sin[gamma])
  /((4/3) r_c q Sin[phi])

Important:
- atildeC and atildeS are derived locally.
- Do not require aC, DeltaC, atildeC0 or atildeS0 as independent parameters.
- Do not add aC0 or DeltaC0 to src/parameters.wl for this implementation.
- chat_plus remains a local helper expression:

  chat_plus = atildeC q Cos[delta_c] + atildeS q Sin[delta_c]

- Because atildeC and atildeS contain q and phi through the denominators, their dependence on q and phi must be handled carefully in the purple constraint.

Project variables:
- R -> use the experimental ratio R from project inputs, or compute it from branching ratios and lifetimes if that is the existing project convention.
- A_CP(pi- K+) -> AcppimKp0, AcppimKpErr
- r -> r0, rErr
- delta -> delta0, deltaErr
- r_c -> rc0, rcErr
- gamma -> gamma0, gammaErr
- q -> qSym when solving symbolically, not q0
- phi -> phiScan or phiNP locally
- rho_c -> rhoc0, rhocErr
- theta_c -> thetac0, thetacErr

If rho_c and theta_c are set to zero or neglected in the working setup, make this explicit with an option rather than silently dropping them.

---

## Avoid self-referential q = f(q)

Do not implement the published expression directly as

qPurple[phi_] := (-quadBc + Sqrt[quadBc^2 - 4 quadAc quadDc])/(2 quadAc)

because quadBc and quadDc contain chat_plus, and chat_plus is proportional to q.

Also, quadDc contains explicit q and q^2 terms.

Safe procedure:

1. Define the full equation symbolically.
2. Substitute

   chat_plus -> atildeC qSym Cos[deltac] + atildeS qSym Sin[deltac]

3. Collect the full equation in qSym.
4. Solve the resulting equation for qSym.
5. Evaluate the solutions numerically as functions of phi.
6. Keep only real positive solutions in the plotting window.

Use a local symbolic variable such as qSym for solving.

Do not use the global q0 as the symbolic solve variable.

---

## Branch diagnostics

Before drawing the final purple band, create diagnostics.

For each phi00 branch, print:

- phi00 value in degrees
- number of real positive q solutions
- phi range where valid points exist
- q range
- whether it contributes to the visible reference region

The reference plot appears to show one dominant visible purple region, but this should emerge from branch filtering rather than being hard-coded.

---

## Suggested implementation structure

Add a new helper file:

src/q_phi_scp_constraint.wl

This file may contain reusable functions such as:

phi00SolutionsFromSCP[scp_, acp_, phid_]

buildSCPPurpleBranchData[...]

solveQForPhi00Branch[...]

Then add a diagnostic script:

scripts/diagnose_q_phi_scp_constraint.wl

This script should:

- load src/init.wl
- load src/q_phi_scp_constraint.wl
- compute phi00 branches
- scan phi over a coarse range
- print branch diagnostics
- optionally export a diagnostic plot

Do not modify the blue/green q-phi implementation in this step.