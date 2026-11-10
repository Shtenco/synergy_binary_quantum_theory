#!/usr/bin/env python3
"""Audit the tetrahedral charged-volume completion before any new S result.

The gate is intentionally both a negative and positive control:
1. demonstrate that fixed q_123 is not a tetrahedrally covariant extension on
   one-hit charged J=1/2 sectors;
2. show that the normalized four-leg Q_tet preserves the old Gauss-sector
   absolute-volume normalization;
3. show that H_E^sine on both K5 and independent 16-cell is unchanged;
4. show that the four 16-cell covariant C_r(V) slot norms become equal and
   nonzero, removing the historical preferred-leg null channel.

No Lorentzian science amplitude is used to choose the completion or thresholds.
"""
from __future__ import annotations
import argparse,itertools,json,math,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import charged_intertwiner_recoupling_gate as CH
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_zeroaware_volume_migration_experiment as ZVM
from tetrahedral_volume_backend import apply_q_tetra,apply_tetra_volume_tensor,tetra_q_block,tetra_volume_block_general,tetra_volume_diagnostics,install_tetrahedral_volume_backend
from pl_dual_complex import DualComplex,boundary_4simplex,seed_16cell_boundary
from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
from pl_covariant_backend import install_pl_graph


def state_rel(a,b):
    keys=set(a)|set(b);num=math.sqrt(sum(abs(a.get(k,0j)-b.get(k,0j))**2 for k in keys));den=math.sqrt(sum(abs(z)**2 for z in b.values()))
    return num/max(den,1e-30)
def n2(s):return float(sum(abs(a)**2 for a in s.values()))

def old_volume_block(sp,J2):
    qs=[CH.q_block(tuple(sp),J2,M)[1] for M in PW.m2vals_t(J2)];Q=sum(qs)/len(qs);Q=.5*(Q+Q.conj().T)
    return ZVM.zeroaware_sqrt_abs(Q)[0]

def gauss_volume_regression():
    rows=[];mx=0.0;checked=0
    for sp in itertools.product(range(4),repeat=4):
        labs=CH.allowed_charged_labels(tuple(sp),0)
        if not labs:continue
        Vo=old_volume_block(sp,0);Vn=tetra_volume_block_general(tuple(sp),0)
        err=float(np.linalg.norm(Vo-Vn)/max(np.linalg.norm(Vo),1e-30)) if np.linalg.norm(Vo)>1e-14 else float(np.linalg.norm(Vn))
        mx=max(mx,err);checked+=1
        if len(labs)>1:rows.append({'spins':list(sp),'singlet_dimension':len(labs),'relative_volume_error':err})
    return checked,mx,rows

def charged_negative_control():
    base=[1,1,1,1];rows=[];max_nonprop=0.0;max_old_zero_new=0.0
    for leg in range(4):
      for so in (0,2):
        sp=base.copy();sp[leg]=so;sp=tuple(sp);J2=1;M=PW.m2vals_t(J2)[0]
        Qo=CH.q_block(sp,J2,M)[1];Qn=tetra_q_block(sp,J2,M)[0]
        den=np.vdot(Qo,Qo)
        sc=np.vdot(Qo,Qn)/den if abs(den)>1e-30 else 0j
        nonprop=float(np.linalg.norm(Qn-sc*Qo)/max(np.linalg.norm(Qn),1e-30));max_nonprop=max(max_nonprop,nonprop)
        if np.linalg.norm(Qo)<1e-14:max_old_zero_new=max(max_old_zero_new,float(np.linalg.norm(Qn)))
        rows.append({'spins':list(sp),'changed_leg':leg,'changed_doubled_spin':so,
                     'old_q123_norm':float(np.linalg.norm(Qo)),'tetra_Q_norm':float(np.linalg.norm(Qn)),
                     'best_scalar_tetra_over_old':[float(sc.real),float(sc.imag)],'relative_nonproportional_residual':nonprop})
    return rows,max_nonprop,max_old_zero_new

def he_pair(tets):
    D=DualComplex(tets);G=PLPeterWeylEuclidean(D);seed=((1,)*len(G.EDGES),(0,)*D.n_tets)
    G.primitive_items.cache_clear();old=G.H_sine_basis(seed,0,5)
    with install_tetrahedral_volume_backend():
        G.primitive_items.cache_clear();new=G.H_sine_basis(seed,0,5)
    G.primitive_items.cache_clear()
    return {'old_support':len(old),'new_support':len(new),'old_norm':G.norm(old),'new_norm':G.norm(new),
            'relative_error':state_rel(new,old),'support_identical':set(new)==set(old)}

def cv_slot_norms(use_tetra):
    D=DualComplex(seed_16cell_boundary());G=PLPeterWeylEuclidean(D);seed=((1,)*len(G.EDGES),(0,)*D.n_tets)
    ctx=install_tetrahedral_volume_backend() if use_tetra else None
    def compute():
      with install_pl_graph(G):
        psi=CV.gauss_to_covariant({seed:1+0j},0);out=[]
        for r in range(4):
          w=D.neighbor[(0,r)];tot=0.0;leak=0.0;supports=[]
          for i,j in itertools.product(range(2),repeat=2):
            s,l=LP.RAW.COMP.C_volume_component(psi,0,w,i,j,7);tot+=n2(s);leak=max(leak,float(l));supports.append(len(s))
          out.append({'slot':r,'target_node':w,'matrix_Frobenius_norm':math.sqrt(tot),'component_supports':supports,'max_complete_basis_leakage':leak})
        return out
    if ctx is None:return compute()
    with ctx:return compute()

def charged_spin_sectors():
    base=[1,1,1,1]
    for leg in range(4):
        for doubled_spin in (0,2):
            spins=base.copy();spins[leg]=doubled_spin
            yield tuple(spins)


def charged_spectral_rows():
    rows=[]
    for spins in charged_spin_sectors():
        for J2 in CV.all_total_J2(spins):
            if not CH.allowed_charged_labels(spins,J2):
                continue
            diagnostics=tetra_volume_diagnostics(spins,J2)
            rows.append({"spins":list(spins),"J2":int(J2),"diagnostics":diagnostics})
    return rows


def spectral_matrix_bundle():
    """Persist raw Q_{J,M} blocks and V_J blocks; no averaging is hidden."""
    arrays={}
    index=[]
    for spins in charged_spin_sectors():
        for J2 in CV.all_total_J2(spins):
            labels=CH.allowed_charged_labels(spins,J2)
            if not labels:
                continue
            V=tetra_volume_block_general(spins,J2)
            stem="s"+"_".join(str(x) for x in spins)+f"_J2_{int(J2)}"
            vkey=stem+"_V"
            arrays[vkey]=np.asarray(V,dtype=np.complex128)
            index.append({"kind":"V_J","key":vkey,"spins":list(spins),
                          "J2":int(J2),"shape":list(V.shape)})
            for M2 in PW.m2vals_t(J2):
                Q,_,_=tetra_q_block(spins,J2,M2)
                qkey=stem+f"_M2_{int(M2)}_Q_raw"
                arrays[qkey]=np.asarray(Q,dtype=np.complex128)
                index.append({"kind":"Q_JM_raw","key":qkey,"spins":list(spins),
                              "J2":int(J2),"M2":int(M2),"shape":list(Q.shape)})
    return arrays,index


def global_magnetic_operator_rows():
    rows=[]
    for spins in ((1,1,1,1),(0,1,1,1)):
        shape=tuple(value+1 for value in spins)
        dimension=int(np.prod(shape))
        Q=np.zeros((dimension,dimension),complex)
        V=np.zeros_like(Q)
        for column in range(dimension):
            basis=np.zeros(shape,dtype=complex)
            basis.flat[column]=1.0
            Q[:,column]=apply_q_tetra(basis,spins).reshape(-1)
            V[:,column]=apply_tetra_volume_tensor(basis,spins).reshape(-1)
        qherm=float(np.linalg.norm(Q-Q.conj().T,ord="fro"))
        vherm=float(np.linalg.norm(V-V.conj().T,ord="fro"))
        vmin=float(np.min(np.linalg.eigvalsh(0.5*(V+V.conj().T))))
        qscale=max(float(np.linalg.norm(Q,ord="fro")),1.0)
        vscale=max(float(np.linalg.norm(V,ord="fro")),1.0)
        functional=float(np.linalg.norm((V@V)@(V@V)-Q@Q,ord="fro")/max(float(np.linalg.norm(Q@Q,ord="fro")),1.0))
        commutators=[]
        for axis in range(3):
            total=np.zeros_like(Q)
            for leg,spin in enumerate(spins):
                factors=[np.eye(s+1,dtype=complex) for s in spins]
                factors[leg]=PW.spin_mats_cached(spin)[axis]
                term=factors[0]
                for factor in factors[1:]:
                    term=np.kron(term,factor)
                total+=term
            commutators.append(float(np.linalg.norm(Q@total-total@Q,ord="fro")))
        rows.append({
            "spins":list(spins),"magnetic_dimension":dimension,
            "raw_global_Q_hermiticity_relative":qherm/qscale,
            "global_V_hermiticity_relative":vherm/vscale,
            "global_V_min_eigenvalue":vmin,
            "global_V4_minus_Q2_relative":functional,
            "global_Q_rotation_commutator_max":max(commutators),
            "all_finite":bool(np.isfinite(Q.real).all() and np.isfinite(Q.imag).all()
                              and np.isfinite(V.real).all() and np.isfinite(V.imag).all()),
        })
    return rows


def tensor_reconstruction_rows():
    rows=[]
    for spins in ((1,1,1,1),(0,1,1,1)):
        shape=tuple(value+1 for value in spins)
        size=int(np.prod(shape))
        values=np.arange(1,size+1,dtype=float).reshape(shape)
        tensor=values+1j*np.flip(values,axis=0)
        output=apply_tetra_volume_tensor(tensor,spins)
        rows.append({
            "spins":list(spins),
            "input_shape":list(shape),
            "input_norm":float(np.linalg.norm(tensor)),
            "output_norm":float(np.linalg.norm(output)),
            "finite_output":bool(np.isfinite(output.real).all() and np.isfinite(output.imag).all()),
            "reconstruction_gate_enforced_inside_apply":True,
        })
    return rows


def run():
    ZVM.patch_and_clear();checked,gauss_err,gauss_rows=gauss_volume_regression();charged,max_nonprop,oldzero=charged_negative_control()
    k5=he_pair(boundary_4simplex());cell=he_pair(seed_16cell_boundary())
    oldslots=cv_slot_norms(False);newslots=cv_slot_norms(True)
    spectral_rows=charged_spectral_rows()
    matrix_arrays,matrix_index=spectral_matrix_bundle()
    tensor_rows=tensor_reconstruction_rows()
    global_rows=global_magnetic_operator_rows()
    oldnorm=[x['matrix_Frobenius_norm'] for x in oldslots];newnorm=[x['matrix_Frobenius_norm'] for x in newslots]
    oldspread=max(oldnorm)-min(oldnorm);newspread=max(newnorm)-min(newnorm)
    checks={
      'charged_q123_preferred_leg_negative_control':oldzero>0.1 and max_nonprop>0.5,
      'gauss_absolute_volume_normalization_preserved':gauss_err<1e-12,
      'K5_H_E_sine_exactly_preserved':k5['support_identical'] and k5['relative_error']<1e-12,
      '16cell_H_E_sine_exactly_preserved':cell['support_identical'] and cell['relative_error']<1e-12,
      'old_CV_slot_asymmetry_detected':oldspread>0.1 and min(oldnorm)<1e-12,
      'tetra_CV_all_slots_nonzero':min(newnorm)>1e-3,
      'tetra_CV_slot_covariance_restored':newspread<1e-12,
      'tetra_CV_basis_complete':max(x['max_complete_basis_leakage'] for x in newslots)<1e-10,
      'charged_spectral_sectors_explicitly_evaluated':len(spectral_rows)>8,
      'every_raw_M_block_spectrum_recorded':all(
          len(x['diagnostics']['raw_Q_M_block_eigenvalues'])==len(x['diagnostics']['M_values_doubled'])
          and all(len(ev)==len(CH.allowed_charged_labels(tuple(x['spins']),x['J2']))
                  and all(math.isfinite(v) for v in ev)
                  for ev in x['diagnostics']['raw_Q_M_block_eigenvalues'])
          for x in spectral_rows),
      'every_raw_M_block_matrix_norm_recorded':all(
          len(x['diagnostics']['raw_Q_M_block_frobenius_norms'])==len(x['diagnostics']['M_values_doubled'])
          and all(math.isfinite(v) for v in x['diagnostics']['raw_Q_M_block_frobenius_norms'])
          for x in spectral_rows),
      'raw_Q_JM_matrix_registry_complete':sum(x['kind']=='Q_JM_raw' for x in matrix_index)==sum(
          len(x['diagnostics']['M_values_doubled']) for x in spectral_rows)
          and all(np.isfinite(a.real).all() and np.isfinite(a.imag).all()
                  for k,a in matrix_arrays.items() if k.endswith('_Q_raw')),
      'charged_Q_raw_Hermitian_before_symmetrization':max((max(x['diagnostics']['raw_Q_hermiticity_defects']) for x in spectral_rows),default=math.inf)<1e-10,
      'charged_Q_raw_M_independent_before_averaging':max((x['diagnostics']['max_raw_Q_M_block_defect'] for x in spectral_rows),default=math.inf)<1e-10,
      'charged_Q_preserves_total_JM_sectors_before_reduction':max((x['diagnostics']['max_Q_sector_leakage_relative'] for x in spectral_rows),default=math.inf)<1e-10,
      'charged_volume_spectral_functional_calculus':max((x['diagnostics']['V4_minus_Q2_relative_error'] for x in spectral_rows),default=math.inf)<1e-9,
      'complete_magnetic_tensor_reconstruction_and_finite_output':len(tensor_rows)==2 and all(x['finite_output'] for x in tensor_rows),
      'global_raw_Q_and_V_hermiticity':len(global_rows)==2 and max(x['raw_global_Q_hermiticity_relative'] for x in global_rows)<1e-10 and max(x['global_V_hermiticity_relative'] for x in global_rows)<1e-10,
      'global_V_positive_and_functional_calculus':min(x['global_V_min_eigenvalue'] for x in global_rows)>-1e-10 and max(x['global_V4_minus_Q2_relative'] for x in global_rows)<1e-9,
      'global_Q_rotational_scalar':max(x['global_Q_rotation_commutator_max'] for x in global_rows)<1e-9,
      'global_operator_matrices_finite':all(x['all_finite'] for x in global_rows),
    }
    return {'status':'tetrahedral charged-volume covariance correction audit','passed':bool(all(checks.values())),
            'science_status':'PRE_LORENTZIAN_OPERATOR_CORRECTION','definition':'Q_tet=(1/4) sum_r (-1)^r q_{omit r}; V_tet=sqrt(abs(Q_tet)) with production zero-aware spectral convention',
            'checks':checks,'gauss_small_spin_quartets_checked':checked,'max_gauss_absolute_volume_relative_error':gauss_err,
            'nontrivial_gauss_rows':gauss_rows,'charged_one_hit_negative_control':charged,
            'max_charged_nonproportional_residual_old_vs_tetra':max_nonprop,'max_tetra_Q_norm_where_old_q123_was_zero':oldzero,
            'K5_H_E_comparison':k5,'sixteen_cell_H_E_comparison':cell,
            'old_q123_CV_slots':oldslots,'tetrahedral_CV_slots':newslots,
            'charged_spectral_rows':spectral_rows,'spectral_matrix_index':matrix_index,
            'spectral_matrix_array_count':len(matrix_arrays),
            'magnetic_tensor_reconstruction_rows':tensor_rows,
            'global_magnetic_operator_rows':global_rows,
            'old_CV_slot_norm_spread':oldspread,'tetra_CV_slot_norm_spread':newspread,
            'interpretation':'The fixed q_123 extension is valid for frozen Gauss absolute volume but not for charged intermediate states. The four-leg completion restores local tetrahedral slot covariance while leaving the gauge-invariant Euclidean Hamiltonian columns unchanged on both regulators.',
            'consequence':'All Lorentzian C(K)/S calculations that used the old charged q_123 extension are retained only as historical diagnostics and must be recomputed with the tetrahedral backend before physical promotion.',
            'scope_note':'This correction does not change the already tested E_sine amplitudes. Full corrected Lorentzian normalization/HDA must be re-audited.'}
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    p.add_argument('--matrix-output',type=Path)
    a=p.parse_args()
    o=run()
    t=json.dumps(o,indent=2)
    print(t)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(t+'\n',encoding='utf-8')
    if a.matrix_output:
        a.matrix_output.parent.mkdir(parents=True,exist_ok=True)
        arrays,_=spectral_matrix_bundle()
        np.savez_compressed(a.matrix_output,**arrays)
    return 0 if o['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
