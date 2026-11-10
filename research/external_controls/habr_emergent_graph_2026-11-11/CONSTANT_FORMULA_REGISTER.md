# Constant and particle formula register — verbatim article transcription

Source: Habr article 1046730, section “Формулы эмерджентных физических констант”, article lines 680–796.  
Rule: each row is an independent audit target. Formula text is transcribed from the article table, not silently replaced with the Python implementation. Numerical verification and dimensional verification are separate columns to be completed.

Common symbols: (K=6), (N) is the claimed graph size. Unless otherwise stated, (ln N) and (ln K) are dimensionless.

| ID | Quantity | Formula as printed | Claimed unit | Audit status |
|---|---|---|---|---|
| P-01 | Planck length (l_P) | (4(ln N)^2ln K/N^{1/3}) | m | TRANSCRIBED |
| P-02 | Planck time (t_P) | (4K^2(ln K)^2/[pi N^{1/3}(ln N)^2]) | s | TRANSCRIBED |
| P-03 | Planck mass (m_P) | (K/[4pi(ln N)^3]) | kg | TRANSCRIBED |
| P-04 | Planck energy (E_P) | (pi(ln N)^5/[4K^3(ln K)^2]) | J | TRANSCRIBED |
| P-05 | Planck temperature (T_P) | (8pi N^{1/3}/(ln N)^4) | K | TRANSCRIBED |
| F-01 | Reduced Planck constant (hbar) | ((ln N)^3/[K N^{1/3}]) | J·s | TRANSCRIBED |
| F-02 | Planck constant (h) | (2pi(ln N)^3/[K N^{1/3}]) | J·s | TRANSCRIBED |
| F-03 | Speed of light (c) | (pi(ln N)^4/[K^2ln K]) | m/s | TRANSCRIBED |
| F-04 | Newton constant (G) | (16pi^3(ln N)^{13}/[K^5ln K,N^{1/3}]) | m³/(kg·s²) | TRANSCRIBED |
| F-05 | Boltzmann constant (k_B) | ((ln N)^8/[8pi^2N^{1/3}]) | J/K | TRANSCRIBED |
| F-06 | Fine-structure constant (alpha) | (2(ln K)^2/[piln N]) | dimensionless | TRANSCRIBED |
| F-07 | Elementary charge (q_e) | (1/[pi K^{3/2}(ln N)^7]) | C | TRANSCRIBED |
| EM-01 | Vacuum permittivity (epsilon_0) | (N^{1/3}/[8pi^3ln K(ln N)^{20}]) | F/m | TRANSCRIBED |
| EM-02 | Vacuum permeability (mu_0) | (8pi K^4(ln K)^3(ln N)^{12}/N^{1/3}) | N/A² | TRANSCRIBED |
| EM-03 | Vacuum impedance (Z_0) | (8K^2pi^2(ln K)^2(ln N)^{16}/N^{1/3}) | Ω | TRANSCRIBED |
| EM-04 | Flux quantum (Phi_0) | (pi^2K^{1/2}(ln N)^{10}/N^{1/3}) | Wb | TRANSCRIBED |
| AT-01 | Rydberg constant (R_infty) | (4(ln N)^3(ln K)^3/[pi K^{3/2}]) | m⁻¹ | TRANSCRIBED |
| AT-02 | Bohr radius (a_0) | (K^{3/2}/[8pi(ln N)^4ln K]) | m | TRANSCRIBED |
| AT-03 | Electron Compton wavelength (lambda_e) | (K^{3/2}ln K/[2pi(ln N)^5]) | m | TRANSCRIBED |
| AT-04 | Proton Compton wavelength (lambda_p) | (2K^{5/2}ln K/[sqrtpi(ln N)^7]) | m | TRANSCRIBED |
| AT-05 | Classical electron radius (r_e) | (K^{3/2}(ln K)^3/[2pi^3(ln N)^6]) | m | TRANSCRIBED |
| L-01 | Electron mass (m_e) | (4pi(ln N)^4/[K^{1/2}N^{1/3}]) | kg | TRANSCRIBED |
| L-02 | Muon mass (m_mu) | (4pi^2(ln N)^5/[Ksqrt3,N^{1/3}]) | kg | TRANSCRIBED |
| L-03 | Tau mass (m_	au) | (sqrtpi K^2(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| Q-01 | Up-quark mass (m_u) | (sqrt3(ln N)^5/[4pi^2N^{1/3}]) | kg | TRANSCRIBED |
| Q-02 | Down-quark mass (m_d) | ((ln N)^5/[Ksqrt3,N^{1/3}]) | kg | TRANSCRIBED |
| Q-03 | Strange-quark mass (m_s) | (K(ln N)^5/[sqrtpi,N^{1/3}]) | kg | TRANSCRIBED |
| Q-04 | Charm-quark mass (m_c) | (2pi^2(ln N)^6/[K^3N^{1/3}]) | kg | TRANSCRIBED |
| Q-05 | Bottom-quark mass (m_b) | (pi(ln N)^6/[Ksqrt3,N^{1/3}]) | kg | TRANSCRIBED |
| Q-06 | Top-quark mass (m_t) | (K^3(ln N)^6/[pi^2N^{1/3}]) | kg | TRANSCRIBED |
| ME-01 | Neutral pion (m_{pi^0}) | (2pi K^3(ln N)^4/N^{1/3}) | kg | TRANSCRIBED |
| ME-02 | Charged pion (m_{pi^pm}) | ((ln N)^6/[4pi^2sqrt2,N^{1/3}]) | kg | TRANSCRIBED |
| ME-03 | Neutral kaon (m_{K^0}) | (sqrt{2pi}(ln N)^6/[4pi^2N^{1/3}]) | kg | TRANSCRIBED |
| ME-04 | Eta meson (m_eta) | (2pi^2(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| ME-05 | Eta-prime (m_{eta'}) | (K^3(ln N)^5/[2pi N^{1/3}]) | kg | TRANSCRIBED |
| ME-06 | Rho meson (m_ho) | (sqrt3pi^{5/2}(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| ME-07 | Omega meson (m_omega) | (2sqrt2pi^2(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| ME-08 | Phi meson (m_phi) | (sqrt{2pi}K^{3/2}(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| ME-09 | K-star (m_{K^*}) | (2pi^{-5/2}(ln N)^6/N^{1/3}) | kg | TRANSCRIBED |
| ME-10 | D-zero (m_{D^0}) | (sqrt{2pi}(ln N)^6/[Ksqrt3,N^{1/3}]) | kg | TRANSCRIBED |
| ME-11 | J/psi (m_{J/psi}) | (8sqrt2pi^2(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| ME-12 | Eta-c (m_{eta_c}) | (4cdot3^3(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| ME-13 | Upsilon(1S) (m_Upsilon) | (sqrt3(ln N)^6/[sqrt2,N^{1/3}]) | kg | TRANSCRIBED |
| ME-14 | B-zero (m_{B^0}) | (pi^2(ln N)^6/[2^{3/2}3^{3/2}N^{1/3}]) | kg | TRANSCRIBED |
| ME-15 | B-c (m_{B_c}) | (4sqrt2cdot9pi(ln N)^5/N^{1/3}) | kg | TRANSCRIBED |
| BA-01 | Proton (m_p) | (sqrtpi(ln N)^6/[K^{3/2}N^{1/3}]) | kg | TRANSCRIBED |
| BA-02 | Lambda (m_Lambda) | (sqrt2(ln N)^6/[pi^2N^{1/3}]) | kg | TRANSCRIBED |
| BA-03 | Sigma-plus (m_{Sigma^+}) | (K(ln N)^6/[4pi^2N^{1/3}]) | kg | TRANSCRIBED |
| BA-04 | Sigma-minus (m_{Sigma^-}) | (K(ln N)^6/[4pi^2N^{1/3}]) | kg | TRANSCRIBED |
| BA-05 | Xi-zero (m_{Xi^0}) | (sqrt{2pi}(ln N)^6/[K^{3/2}N^{1/3}]) | kg | TRANSCRIBED |
| BA-06 | Omega-minus (m_{Omega^-}) | (pi(ln N)^6/[K^{3/2}N^{1/3}]) | kg | TRANSCRIBED |
| BA-07 | Lambda-c-plus (m_{Lambda_c}) | (sqrtpi(ln N)^6/[KN^{1/3}]) | kg | TRANSCRIBED |
| BA-08 | Xi-c-plus (m_{Xi_c}) | ((ln N)^6/[pi N^{1/3}]) | kg | TRANSCRIBED |
| BA-09 | Omega-c-zero (m_{Omega_c}) | (K(ln N)^6/[pi^{5/2}N^{1/3}]) | kg | TRANSCRIBED |
| BA-10 | Lambda-b-zero (m_{Lambda_b}) | (sqrtpi(ln N)^6/[sqrt K,N^{1/3}]) | kg | TRANSCRIBED |
| BO-01 | W boson (m_W) | (2pi^3(ln N)^6/[KN^{1/3}]) | kg | TRANSCRIBED |
| BO-02 | Z boson (m_Z) | (4pi^{5/2}(ln N)^6/[KN^{1/3}]) | kg | TRANSCRIBED |
| BO-03 | Higgs boson (m_H) | (4pi^2(ln N)^6/[K^{1/2}N^{1/3}]) | kg | TRANSCRIBED |
| BO-04 | Higgs vacuum scale (v_H) | (8pi^{3/2}(ln N)^6/[sqrt2,N^{1/3}]) | kg (as printed) | TRANSCRIBED |
| CO-01 | Cosmological constant (Lambda) | ((ln N)^{12}/[sqrtpi,N^{2/3}]) | m⁻² | TRANSCRIBED |
| CO-02 | Einstein constant (kappa) | (128K^3(ln K)^3/[(ln N)^3N^{1/3}]) | s²/(kg·m) | TRANSCRIBED |
| CO-03 | Connection energy (E_{conn}) | (8pi(ln N)^{10}(ln K)^2/[K^{9/2}N^{1/3}]) | J | TRANSCRIBED |
| TL-01 | Neutron lifetime (	au_n) | (piln N) | s | TRANSCRIBED |
| TL-02 | Muon lifetime (	au_mu) | (ln K/[Ksqrt3(ln N)^2]) | s | TRANSCRIBED |
| TL-03 | Tau lifetime (	au_	au) | (1/[2(ln N)^5]) | s | TRANSCRIBED |
| TL-04 | Charged-pion lifetime (	au_{pi^pm}) | (K^2sqrt2pi/(ln N)^4) | s | TRANSCRIBED |
| TL-05 | Charged-kaon lifetime (	au_{K^pm}) | (4/[K^{3/2}(ln N)^3]) | s | TRANSCRIBED |
| TL-06 | D-plus lifetime (	au_{D^pm}) | (1/[sqrtpi K^{5/2}(ln N)^4]) | s | TRANSCRIBED |
| TL-07 | B-plus lifetime (	au_{B^pm}) | (piln K/[2(ln N)^5]) | s | TRANSCRIBED |
| TL-08 | Lambda-b lifetime (	au_{Lambda_b}) | (sqrt2ln K/(ln N)^5) | s | TRANSCRIBED |
| TL-09 | D-zero lifetime (	au_{D^0}) | (ln K/[2pi^2K^2(ln N)^4]) | s | TRANSCRIBED |
| R-01 | Proton/electron mass ratio | ((ln N)^2/[4sqrtpi K]) | dimensionless | TRANSCRIBED |
| R-02 | Muon/electron mass ratio | (pisqrt Kln N/sqrt3) | dimensionless | TRANSCRIBED |
| R-03 | Tau/electron mass ratio | (K^{5/2}ln N/[4sqrtpi]) | dimensionless | TRANSCRIBED |
| R-04 | W/Z mass ratio | (sqrtpi/2) | dimensionless | TRANSCRIBED |
| R-05 | Higgs/W mass ratio | (2sqrt K/pi) | dimensionless | TRANSCRIBED |
| R-06 | Planck/electron mass ratio | (K^{3/2}N^{1/3}/[16pi^2(ln N)^7]) | dimensionless | TRANSCRIBED |

## Mandatory follow-up

- Recompute every row using (a) article's own (N_{math}), (b) article's (N_{phys}), (c) code's hard-coded JSON (N=4.1977\times10^{121}), and (d) values returned by the code's own fitting procedure. Do not mix these cases.
- Compare every row to the actual current metrology source and preserve uncertainty, renormalization scheme, and scale. In particular, quark masses are scheme/scale dependent; unstable resonance masses and lifetimes need correct definitions.
- Verify every claimed unit. As written, all RHS expressions are dimensionless, so a dimensional scale is missing unless an explicit natural-unit convention and conversion map are supplied.
- Cross-check each article row against the Python implementation. Any discrepancy gets a separate SOURCE-MISMATCH record.
