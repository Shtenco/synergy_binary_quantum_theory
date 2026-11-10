#!/usr/bin/env python3
"""Independent transcription of physics_constants/stress_test.py.

No SciPy, no fitting and no network access. Evaluates the published formula
dictionary at four explicitly separated N conventions. Source blob SHA:
b617ac8f42bec5e7916dfefce576665bb68e81fc
"""
from __future__ import annotations
import math

K = 6.0
PI = math.pi
LNK = math.log(K)

TARGET = {
    "ħ":1.054571817e-34, "h":6.62607015e-34,
    "t_P":5.391247e-44, "l_P":1.616255e-35, "m_P":2.176434e-8,
    "E_P":1.956082e9, "T_P":1.416784e32,
    "c":299792458, "G":6.67430e-11, "k_B":1.380649e-23,
    "α":1/137.035999084, "m_e":9.1093837015e-31, "r_e":2.8179402853e-15,
    "ep_0":8.8725366415e-12, "mu_0":1.25663706127e-6, "q_e":1.602176634e-19,
    "m_proton":1.67262192e-27, "m_muon":1.883531627e-28, "m_tau":3.167e-27,
    "m_pi_meson":2.4880888e-28, "m_pi0_meson":2.40609e-28,
    "m_k0_meson":8.801929e-28, "m_DT":3.3435837724e-27,
    "m_Λ_barion":1.9901611e-27, "m_z_bozon":1.62614e-25,
    "m_w_bozon":1.43362e-25, "m_Higgs":2.23319e-25,
    "m_D0":3.32479e-27, "m_J_ψ":5.52061e-27, "m_eta":9.767732e-28,
    "m_Υ_1S":1.68715e-26, "m_qu_u":2.1650e-30, "m_qu_d":4.7915e-30,
    "m_qu_s":1.66e-28, "m_qu_c":1.27e-27, "m_qu_b":4.180e-27,
    "m_qu_t":3.04e-25, "RIDBERG":1.097373e7, "bor_radius":5.29177210903e-11,
    "impedance":376.730313, "Φ0_magnetic_stream":2.06783366752e-15,
    "m_proton_to_m_electron":1836.152673426, "m_tau_m_electron":3477,
    "m_W_to_m_Z":0.8815, "m_plank_to_m_e":2.389e22,
    "compton_e":2.426e-12, "compton_proton":1.32140985396e-15,
    "m_Higgs_to_m_W":1.558, "Lambda_cosmo":1.08929e-52,
    "Einstein_constant":2.07664746e-43, "vacuum_higgs":4.388471e-25,
    "mu_lifetime":2.1969811e-6, "tau_lifetime":2.903e-13,
    "pion_lifetime":2.6033e-8, "neutron_lifetime":877.8,
    "kaon_lifetime":1.2380e-8, "D_plus_lifetime":1.040e-12,
    "B_plus_lifetime":1.638e-12, "Λ_b_lifetime":1.471e-12,
    "D0_lifetime":4.101e-13, "h2_connection_energy":2.178872e-18,
    "Sigma_plus":2.11933e-27, "m_neitrino":1.783e-36,
    "Ksi_0":2.34532e-27, "omega_minus":2.9859e-27,
    "Ksi_plus":4.3995e-27, "Omega0_c":4.808e-27, "lambda_B0":1.0023e-26,
    "Ksi_minus":2.358e-27, "Sigma_minus":2.132e-27,
    "phi_meson":1.819e-27, "omega_meson":1.394e-27,
    "eta_shtrih":1.7086e-27, "rho_meson":1.49e-27,
    "K_star_meson":1.59e-27, "eta_c_meson":5.319e-27,
    "h_c_meson":6.285e-27, "delta_meson":2.196e-27,
    "B_c_meson":1.1185e-26, "Ksi_pp_b_baryon":6.543e-27,
    "K_plus":8.802e-28, "B_meson":9.413e-27,
}


def formulas(n: float) -> dict[str, float]:
    ln = math.log(n)
    n13 = n ** (1.0 / 3.0)
    return {
        "ħ": ln**3/(K*n13),
        "h": 2*PI*ln**3/(K*n13),
        "t_P": 4*K**2*LNK**2/(PI*n13*ln**2),
        "l_P": 4*ln**2*LNK/n13,
        "m_P": K/(4*PI*ln**3),
        "E_P": ln**5*PI/(4*K**3*LNK**2),
        "T_P": 8*PI*n13/ln**4,
        "c": PI*ln**4/(K**2*LNK),
        "G": 16*PI**3*ln**13/(K**5*LNK*n13),
        "k_B": ln**8/(8*PI**2*n13),
        "α": 2*LNK**2/(PI*ln),
        "m_e": 4*PI*ln**4/(math.sqrt(K)*n13),
        "m_proton": math.sqrt(PI)*ln**6/(K**1.5*n13),
        "ep_0": n13/(8*PI**3*LNK*ln**20),
        "mu_0": 8*PI*K**4*LNK**3*ln**12/n13,
        "q_e": 1/(PI*K**1.5*ln**7),
        "m_muon": 4*PI**2*ln**5/(K*math.sqrt(3)*n13),
        "m_tau": math.sqrt(PI)*ln**5*K**2/n13,
        "m_pi_meson": ln**6/(4*PI**2)*n**(-1/3)/math.sqrt(2),
        "m_pi0_meson": 2*PI*K**3*ln**4/n13,
        "m_k0_meson": ln**6/(4*PI**2)*math.sqrt(2*PI)/n13,
        "m_DT": ln**6*math.sqrt(2*PI)/(K*math.sqrt(3)*n13),
        "m_Λ_barion": ln**6*n**(-1/3)*math.sqrt(2)/PI**2,
        "m_z_bozon": ln**6*4*PI**2.5/(n13*K),
        "m_w_bozon": 2*PI**3*ln**6/(n13*K),
        "m_Higgs": ln**6*4*PI**2/(n13*K**0.5),
        "m_D0": ln**6*math.sqrt(2*PI)/(n13*K*math.sqrt(3)),
        "m_J_ψ": ln**5*8*PI**2*math.sqrt(2)/n13,
        "m_eta": ln**5*2*PI**2/n13,
        "m_Υ_1S": ln**6*math.sqrt(3)/(math.sqrt(2)*n13),
        "m_qu_u": ln**5*math.sqrt(3)/(4*PI**2*n13),
        "m_qu_d": ln**5/(K*math.sqrt(3)*n13),
        "m_qu_s": ln**5*K/(math.sqrt(PI)*n13),
        "m_qu_c": ln**6*2*PI**2/(K**3*n13),
        "m_qu_b": ln**6*PI/(K*math.sqrt(3)*n13),
        "m_qu_t": ln**6*K**3/(PI**2*n13),
        "RIDBERG": 4*ln**3*LNK**3/(PI*K**1.5),
        "r_e": LNK**3*K**1.5/(2*PI**3*ln**6),
        "bor_radius": K**1.5/(8*PI*ln**4*LNK),
        "impedance": 8*K**2*PI**2*LNK**2*ln**16/n13,
        "Φ0_magnetic_stream": ln**10*PI**2*K**0.5/n13,
        "m_proton_to_m_electron": ln**2/(4*math.sqrt(PI)*K),
        "m_tau_m_electron": K**2.5*ln/(4*math.sqrt(PI)),
        "m_W_to_m_Z": math.sqrt(PI)/2,
        "m_plank_to_m_e": K**1.5*n13/(16*PI**2*ln**7),
        "compton_e": K**1.5*LNK/(2*PI*ln**5),
        "compton_proton": 2*K**2.5*LNK/(math.sqrt(PI)*ln**7),
        "m_Higgs_to_m_W": 2*K**0.5/PI,
        "Lambda_cosmo": ln**12/(math.sqrt(PI)*n**(2/3)),
        "Einstein_constant": 128*K**3*LNK**3/(ln**3*n13),
        "vacuum_higgs": ln**6*8*PI**1.5/(math.sqrt(2)*n13),
        "mu_lifetime": LNK/(K*math.sqrt(3)*ln**2),
        "tau_lifetime": 1/(2*ln**5),
        "pion_lifetime": K**2*math.sqrt(2)*PI/ln**4,
        "D_plus_lifetime": 1/(math.sqrt(PI)*K**2.5*ln**4),
        "neutron_lifetime": PI*ln,
        "kaon_lifetime": 4/(K**1.5*ln**3),
        "B_plus_lifetime": LNK*PI/(2*ln**5),
        "Λ_b_lifetime": LNK*math.sqrt(2)/ln**5,
        "D0_lifetime": LNK/(2*PI**2*K**2*ln**4),
        "h2_connection_energy": 8*PI*ln**10*LNK**2/(K**4.5*n13),
        "Sigma_plus": K*ln**6/(4*PI**2*n13),
        "m_neitrino": ln**2*n**(-1/3)*math.sqrt(2)/LNK,
        "Ksi_0": ln**6*math.sqrt(2*PI)*n**(-1/3)/K**1.5,
        "omega_minus": ln**6*PI*n**(-1/3)/K**1.5,
        "Ksi_plus": ln**6/(PI*n13),
        "Omega0_c": K**0.5*ln**5*4*PI**2/n13,
        "lambda_B0": ln**6*math.sqrt(PI)/(K**0.5*n13),
        "Ksi_minus": ln**6*math.sqrt(2*PI)*n**(-1/3)/K**1.5,
        "Sigma_minus": ln**6*K*n**(-1/3)/(4*PI**2),
        "phi_meson": ln**5*math.sqrt(2*PI)*K**1.5/n13,
        "omega_meson": ln**5*2*PI**2*math.sqrt(2)/n13,
        "eta_shtrih": ln**5*K**3/(2*PI*n13),
        "rho_meson": math.sqrt(3)*PI**2.5*ln**5/n13,
        "K_star_meson": 2*PI**(-2.5)*ln**6/n13,
        "B_meson": PI**2*ln**6/(2**1.5*3**1.5*n13),
        "eta_c_meson": (math.sqrt(2))**4*(math.sqrt(3))**6*ln**5/n13,
        "h_c_meson": math.sqrt(2)*ln**6/(PI*n13),
        "delta_meson": ln**6/(2*PI*n13),
        "B_c_meson": (math.sqrt(2))**6*(math.sqrt(3))**4*PI*ln**5/n13,
        "B_s_meson": (math.sqrt(2))**4*math.sqrt(3)*ln**6/n13,
        "Ksi_pp_b_baryon": math.sqrt(2)*ln**6/(3*n13),
        "K_plus": ln**6*n**(-1/3)/(2*PI)**1.5,
    }


def main():
    n_math = 4.475947e121
    ln_n_math = math.log(n_math)
    n_script = math.exp(ln_n_math - PI*K/ln_n_math)
    n_geo = math.exp((K-math.log(K))/(1/3-1/PI))
    n_json = 4.1977e121
    n_variation = 4.2e121
    cases = {
        "geometric_exact": n_geo,
        "stress_test_script": n_script,
        "mass_relation_json": n_json,
        "variation_script": n_variation,
    }
    print(f"targets={len(TARGET)}")
    for name, n in cases.items():
        out = formulas(n)
        missing = sorted(set(TARGET)-set(out))
        extra = sorted(set(out)-set(TARGET))
        errors = {k: abs(out[k]/TARGET[k]-1)*100 for k in TARGET if k in out}
        print(f"\n[{name}] N={n:.12e} lnN={math.log(n):.12f}")
        print(f"formula_outputs={len(out)} missing={missing} unmatched_outputs={extra}")
        print(f"mean_relative_error_pct={sum(errors.values())/len(errors):.8f}")
        print(f"median_relative_error_pct={sorted(errors.values())[len(errors)//2]:.8f}")
        for threshold in (0.1, 0.5, 1.0, 5.0):
            print(f"below_{threshold:g}_pct={sum(e<threshold for e in errors.values())}/{len(errors)}")
        for key in sorted(TARGET):
            if key in out:
                print(f"{key}\t{out[key]:.12e}\t{TARGET[key]:.12e}\t{errors[key]:.8f}%")
    print("\nWARNING: These target values are copied from the author's dictionary,")
    print("not independently validated metrology values. Formula evaluation does")
    print("not establish physical derivation or correct units.")


if __name__ == "__main__":
    main()
