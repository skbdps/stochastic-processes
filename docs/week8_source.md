# Week 8 source generation and analytical validation

This module implements the [fixed E4 protocol](week8_exploratory_protocol.md)
and its [mathematical derivation](week8_mathematics.md). Its tests and separate
source-only validation check the sampling mechanism. They do not supply E4
estimator results or change any numerical inference rule.

## Source contract

Every canonical cell is independently addressed by the complete scientific
cell ID and replication number, under namespace 8. The first spawned child
supplies phase and gap randomness; the second supplies stationary OU X0 and
all exact-transition innovations. Execution order, workers and partition size
do not change these streams. A source sample is generated once and supplied
unchanged to exact, PFML and Euler inference.

Time zero is an event. An MMPP sequence starts with the intensity-biased Palm
phase probabilities and retains the actual phase at each arrival. An IID
sequence starts every gap with a new independent Palm phase. Both mechanisms
allow any number of hidden switches before the next arrival. Thus their
entire one-gap distributions agree at each intensity ratio. A mixture of
single exponentials without internal switches would not implement this
control. At ratio one, both gap laws are independent exponential; the regular
reference remains a separate constant-gap design.

The floor is applied after generating raw gaps and latent phase transitions.
It never advances the hidden chain through the added padding and never
renormalizes the gaps. Only ratio 20 triggers the floor in the adopted grid.
The raw vector, the transformed vector, and differences of stored cumulative
times are distinct numerical objects: raw/transformed moments describe their
respective vectors, while `realized_*` diagnostics describe the actual
intervals supplied to the OU simulator. Nonpositive or nonfinite stored
intervals cause a retained source failure; they are not repaired or redrawn.

## Analytical checks

For row-vector Palm law alpha, no-arrival generator D0, A=-D0 and diagonal
arrival matrix Lambda, the raw survival is alpha exp(D0 t) 1 and its kth
moment is k! alpha A^(-k) 1. For floor f, the transformed first and second
moments are the formulas in the mathematical appendix. The first-moment
transition kernel is R_f=fP+exp(D0 f) A^(-2) Lambda, where P=A^(-1)Lambda.
Adjacent MMPP products have expectation alpha R_f R_f 1; independent-control
products have expectation E[Y]^2. The transformed CDF is zero below f and
equals the raw CDF at and above f, including the atom at the floor.

Deterministic tests cover calibration, exact regular gaps, explicit
arrival/switch traces, retained versus reset phase, floor semantics, strict
inputs, immutable samples and deterministic source addresses. The separately
seeded stochastic validation uses exactly 256 independent length-999
sequences in each of six stochastic cells, with no OU paths or fits.
Per-sequence moments and empirical CDFs use 999 gaps, and adjacent products
use 998 pairs. Monte Carlo standard errors use the sample standard deviation
across the 256 independent sequences divided by sqrt(256). They do not treat
dependent gaps as independent observations. The fixed acceptance rule is
absolute discrepancy at most six MCSEs plus 1e-10.

Only the first phase in each independent sequence enters the initial-phase
check. Its reference standard error is binomial with N=256. Clipping is
reported without an empirical rare-event normal gate. Failed validation is
preserved and investigated without changing seeds, replication count or
tolerances to obtain acceptance.

Implementation review, exact commands and observed validation status are
recorded in [the implementation checkpoints](week8_implementation.md).
