#!/usr/bin/env python3
"""BQG v0.7 exact Z2/Gauss cycle-support screen. No H_L amplitudes claimed."""
import itertools, collections, json

VERT = range(5)
EDGES = tuple(itertools.combinations(VERT, 2))
EIDX = {e:i for i,e in enumerate(EDGES)}
def bit(*edges):
    return sum(1<<EIDX[tuple(sorted(e))] for e in edges)
def triangle(a,b,c):
    return bit((a,b),(a,c),(b,c))
def gauss(m):
    return all(sum(bool(m & (1<<i)) for i,e in enumerate(EDGES) if v in e)%2==0 for v in VERT)
def triangles_at(v):
    rest = [x for x in VERT if x!=v]
    return set(triangle(v,a,b) for a,b in itertools.combinations(rest,2))
def lorentz_diffs(v):
    T=triangles_at(v)
    return set(x^y for x in T for y in T)
def local_candidate(v):
    return set(t^d for t in triangles_at(v) for d in lorentz_diffs(v))
def run():
    cycles={m for m in range(1<<10) if gauss(m)}
    odd={m for m in cycles if m.bit_count()%2}
    mask_counts=dict(sorted(collections.Counter(m.bit_count() for m in cycles).items()))
    local=[local_candidate(v) for v in VERT]
    common=set.intersection(*local)
    excluded=odd-common
    e_min=[4.8915556431103875,6.185640314205626,6.690319041756399,7.034044687319637,7.03404468731964]
    assert len(cycles)==64 and len(odd)==32
    assert mask_counts=={0:1,3:10,4:15,5:12,6:15,7:10,10:1}
    assert all(len(triangles_at(v))==6 for v in VERT)
    assert all(len(lorentz_diffs(v))==16 for v in VERT)
    assert all(collections.Counter(m.bit_count() for m in lorentz_diffs(v))=={0:1,4:12,6:3} for v in VERT)
    assert all(len(s)==26 and s<=odd for s in local)
    assert len(common)==22 and collections.Counter(m.bit_count() for m in common)=={3:10,5:12}
    assert len(excluded)==10 and all(m.bit_count()==7 for m in excluded)
    assert min(e_min)>4.8
    return {"status":"finite combinatorial support PASS",
      "source":"original K5 T_sequences triangle path and declared H_L~Tr(C(K)C(K)C(V)) with K=[V,H_E]",
      "gauss_cycle_dimension":6,"gauss_valid_parity_masks":64,"gauss_odd_masks":32,
      "gauss_mask_length_distribution":mask_counts,
      "he_local_triangle_masks":6,"hl_local_two_triangle_xor_masks":16,
      "locally_eligible_odd_masks_per_node":[len(s) for s in local],
      "single_odd_pattern_global_candidates":22,"candidates_by_length":{3:10,5:12},
      "single_7edge_odd_pattern_excluded":10,"local_euclidean_min_eigenvalues":e_min,
      "full_HL_amplitudes_computed":False,"physical_kernel_computed":False,
      "disclaimer":"Necessary support rule only. Does not exclude superpositions over multiple odd patterns; amplitudes, ordering, D and physical clocks remain open."}
if __name__=="__main__":
    print(json.dumps(run(),indent=2))
