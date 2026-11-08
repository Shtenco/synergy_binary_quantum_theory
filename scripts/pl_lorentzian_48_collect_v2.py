#!/usr/bin/env python3
"""Provenance-locked collector for tetrahedral-volume PL Lorentzian V2 terms."""
from __future__ import annotations
import argparse,hashlib,json,re,subprocess,traceback
from pathlib import Path
import numpy as np
import pl_lorentzian_48_collect as BASE
VERSION='tetrahedral-charged-volume-v2'

def _sha(value):
    return isinstance(value,str) and re.fullmatch(r'[0-9a-f]{64}',value) is not None


def _int_list(value):
    return isinstance(value,list) and all(type(item) is int for item in value)


def validate_term_record(path,metadata,expected_source=0):
    """Validate one V2 metadata/NPZ pair before the 48-term reduction."""
    if not isinstance(metadata,dict) or metadata.get('passed') is not True:
        raise RuntimeError('V2_TERM_NOT_PASSED_BOOLEAN_TRUE')
    mode=metadata.get('mode')
    index=metadata.get('index')
    source=metadata.get('source_node')
    coefficient=metadata.get('PL_epsilon_coefficient')
    if mode not in ('forward','adjoint'):
        raise RuntimeError('V2_TERM_BAD_MODE')
    if type(index) is not int or index not in range(24):
        raise RuntimeError('V2_TERM_BAD_INTEGER_INDEX')
    if type(source) is not int or source != expected_source:
        raise RuntimeError('V2_TERM_BAD_SOURCE_NODE')
    if type(coefficient) is not int or coefficient not in (-1,1):
        raise RuntimeError('V2_TERM_BAD_ORIENTATION_COEFFICIENT')
    jmax=metadata.get('Jmax')
    if (not isinstance(jmax,(int,float)) or isinstance(jmax,bool)
            or float(jmax) != 3.5):
        raise RuntimeError('V2_TERM_BAD_CUTOFF')
    checks=metadata.get('checks')
    if (not isinstance(checks,dict) or not checks
            or not all(value is True for value in checks.values())):
        raise RuntimeError('V2_TERM_CHECKS_NOT_BOOLEAN_TRUE')
    npz_path=Path(path).with_suffix('.npz')
    digest=metadata.get('npz_sha256')
    if not _sha(digest) or not npz_path.is_file():
        raise RuntimeError('V2_TERM_NPZ_DIGEST_MISSING')
    if hashlib.sha256(npz_path.read_bytes()).hexdigest() != digest:
        raise RuntimeError('V2_TERM_NPZ_DIGEST_MISMATCH')
    with np.load(npz_path,allow_pickle=False) as data:
        if set(data.files) != {'spins','Ks','amp'}:
            raise RuntimeError('V2_TERM_NPZ_BAD_KEYS')
        spins,Ks,amp=data['spins'],data['Ks'],data['amp']
        if (spins.shape != (len(amp),32) or Ks.shape != (len(amp),16)
                or amp.ndim != 1 or spins.dtype.kind not in 'iu'
                or Ks.dtype.kind not in 'iu' or amp.dtype.kind != 'c'
                or not np.isfinite(amp).all()):
            raise RuntimeError('V2_TERM_NPZ_BAD_SHAPE_DTYPE_OR_AMPLITUDE')
        support=metadata.get('gauss_support')
        if type(support) is not int or support != len(amp):
            raise RuntimeError('V2_TERM_NPZ_SUPPORT_MISMATCH')
        seen=set()
        for spin_row,K_row in zip(spins,Ks):
            key=tuple(map(int,spin_row)),tuple(map(int,K_row))
            if key in seen:
                raise RuntimeError('V2_TERM_NPZ_DUPLICATE_KEY')
            seen.add(key)
    return mode,index,source


def preflight(root):
    rows=[]
    for path in Path(root).rglob('term_*.json'):
        metadata=json.loads(path.read_text(encoding='utf-8'))
        rows.append((path,metadata))
    if len(rows)!=48:
        raise RuntimeError(f'need exactly 48 V2 term metadata files, got {len(rows)}')
    bad=[str(path) for path,metadata in rows
         if metadata.get('operator_version')!=VERSION]
    if bad:
        raise RuntimeError(f'non-V2/superseded term metadata present: {bad[:8]}')

    from pl_dual_complex import DualComplex,seed_16cell_boundary
    import pl_lorentzian_triple_worker as W
    keys=[]
    commits=set()
    backend_ids=set()
    backend_shas=set()
    source_nodes=set()
    D=DualComplex(seed_16cell_boundary())
    for path,metadata in rows:
        mode,index,source=validate_term_record(path,metadata,expected_source=0)
        if path.stem != f'term_{mode}_{index}':
            raise RuntimeError('V2_TERM_METADATA_FILENAME_MISMATCH')
        omit,base,perm,targets,coefficient=W.ordered_spec(D,source,index)
        if metadata.get('PL_epsilon_coefficient') != coefficient:
            raise RuntimeError('V2_TERM_ORIENTATION_COEFFICIENT_MISMATCH')
        if (type(metadata.get('omitted_local_slot')) is not int
                or metadata.get('omitted_local_slot') != omit
                or not _int_list(metadata.get('base_local_slots'))
                or metadata.get('base_local_slots') != list(base)
                or not _int_list(metadata.get('permuted_local_slots'))
                or metadata.get('permuted_local_slots') != list(perm)
                or not _int_list(metadata.get('ordered_target_nodes'))
                or metadata.get('ordered_target_nodes') != list(targets)):
            raise RuntimeError('V2_TERM_ORDERED_PATH_METADATA_MISMATCH')
        keys.append((mode,index))
        source_nodes.add(source)
        commits.add(metadata.get('source_commit'))
        backend_ids.add(metadata.get('volume_backend'))
        backend_shas.add(metadata.get('volume_backend_sha256'))

    expected={(mode,index) for mode in ('forward','adjoint') for index in range(24)}
    if len(set(keys))!=48 or set(keys)!=expected:
        raise RuntimeError('V2_TERM_ORBIT_NOT_EXACTLY_24_BY_2')
    if source_nodes!={0}:
        raise RuntimeError('V2_TERM_SOURCE_NODES_MISSING_OR_MIXED')
    if len(commits)!=1 or not all(
            isinstance(value,str) and re.fullmatch(r'[0-9a-f]{40}',value)
            for value in commits):
        raise RuntimeError('V2_TERM_SOURCE_COMMITS_MISSING_OR_MIXED')
    if backend_ids!={'tetrahedral-volume-v1'} or len(backend_shas)!=1:
        raise RuntimeError('V2_TERM_VOLUME_BACKEND_IDENTITY_MISSING_OR_MIXED')
    if not all(_sha(value) for value in backend_shas):
        raise RuntimeError('V2_TERM_BACKEND_SHA_MALFORMED')
    try:
        actual_head=subprocess.check_output(
            ['git','rev-parse','HEAD'],text=True,timeout=10).strip()
        dirty=subprocess.check_output(
            ['git','status','--porcelain','--untracked-files=no'],text=True,timeout=10)
    except (OSError,subprocess.SubprocessError) as exc:
        raise RuntimeError('V2_COLLECTOR_SOURCE_IDENTITY_NOT_VERIFIABLE') from exc
    if actual_head not in commits:
        raise RuntimeError('V2_COLLECTOR_SOURCE_COMMIT_MISMATCH')
    if dirty.strip():
        raise RuntimeError('V2_COLLECTOR_TRACKED_WORKTREE_NOT_CLEAN')
    backend_file=Path(__file__).with_name('tetrahedral_volume_backend.py')
    actual_backend_sha=hashlib.sha256(backend_file.read_bytes()).hexdigest()
    if backend_shas!={actual_backend_sha}:
        raise RuntimeError('V2_COLLECTOR_BACKEND_SHA_MISMATCH')
    return rows

def run(root):
    preflight(root);L,Ld,S,out=BASE.run(root)
    out['operator_version']=VERSION
    out['charged_volume_definition']='V_tet=sqrt(abs((1/4) sum_r (-1)^r q_hat_r)) with production zero-aware spectrum'
    out['superseded_input_rejected']=True
    out['science_status']='AMPLITUDE_PRECURSOR_S_NODE0_V2'
    return L,Ld,S,out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--json-output',type=Path,required=True);p.add_argument('--state-output',type=Path,required=True);a=p.parse_args()
    try:L,Ld,S,out=run(a.root);code=0 if out['passed'] else 1
    except Exception as exc:L=Ld=S={};out={'status':'V2 collector exception','passed':False,'operator_version':VERSION,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    a.json_output.parent.mkdir(parents=True,exist_ok=True);a.json_output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8');BASE.save_bundle(a.state_output,L,Ld,S);print(json.dumps(out,indent=2,sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
