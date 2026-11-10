#!/usr/bin/env python3
"""Inventory tiny persisted [3,2] master-map key metadata across recovery runs.

This script never touches the depth-6 shell, orbit construction, multiplicity,
Jucys selectors, Hamiltonian, or raw master matrices. It consumes only already
extracted BQG_MIXED_MASTER_KEY_METADATA artifacts, deduplicates shard IDs by
latest artifact creation time, and reports the fail-closed actual-q coverage
needed before any global uniqueness schedule is allowed.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import stat
import sys
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'proof_tools'))
from safe_pickle import safe_load_builtin_pickle

import build_master_key_peeling_schedule as S

NAME_RE = re.compile(r'^bqg-mixed-32-shard-(\d+)-keymeta(?:-targeted-retry)?$')
SALVAGE_RE = re.compile(r'^bqg-depth6-32-stagea-keymeta-salvage-\d+$')
GITHUB_API_HOST = 'api.github.com'


def github_headers(token: str) -> dict[str, str]:
    return {
        'Authorization': f'Bearer {token}',
        'Accept': 'application/vnd.github+json',
        'X-GitHub-Api-Version': '2022-11-28',
        'User-Agent': 'bqg-keymeta-inventory',
    }


def api_json(url: str, token: str) -> dict:
    req = urllib.request.Request(url, headers=github_headers(token))
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


class StripAuthRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Never forward the GitHub bearer token to signed blob-storage redirects."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is None:
            return None
        parsed = urllib.parse.urlparse(newurl)
        if parsed.scheme.lower() != 'https':
            raise urllib.error.HTTPError(
                req.full_url, code, 'refusing non-HTTPS artifact redirect', headers, fp
            )
        old_host = urllib.parse.urlparse(req.full_url).hostname
        new_host = parsed.hostname
        if old_host == GITHUB_API_HOST and new_host != GITHUB_API_HOST:
            new.remove_header('Authorization')
        return new


def download(url: str, token: str, out: Path) -> tuple[int, str]:
    req = urllib.request.Request(url, headers=github_headers(token))
    opener = urllib.request.build_opener(StripAuthRedirectHandler())
    digest = hashlib.sha256()
    size = 0
    with opener.open(req, timeout=180) as response, out.open('xb') as f:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
            digest.update(chunk)
            size += len(chunk)
    return size, digest.hexdigest()


def list_run_artifacts(repo: str, run_id: int, token: str) -> list[dict]:
    all_items = []
    page = 1
    while True:
        x = api_json(
            f'https://api.github.com/repos/{repo}/actions/runs/{run_id}/artifacts?per_page=100&page={page}',
            token,
        )
        items = x.get('artifacts', [])
        all_items.extend(items)
        if len(items) < 100:
            break
        page += 1
    return all_items


def is_candidate_keymeta_artifact_name(name: str) -> bool:
    return bool(NAME_RE.match(name) or SALVAGE_RE.match(name))


def choose_latest_keymeta(repo: str, run_ids: list[int], token: str) -> dict[int, dict]:
    chosen: dict[int, dict] = {}
    for run_id in run_ids:
        for a in list_run_artifacts(repo, run_id, token):
            m = NAME_RE.match(str(a.get('name', '')))
            if not m or a.get('expired'):
                continue
            sid = int(m.group(1))
            rec = dict(a)
            rec['source_keymeta_run_id'] = run_id
            old = chosen.get(sid)
            if old is None or str(rec.get('created_at', '')) > str(old.get('created_at', '')):
                chosen[sid] = rec
    return chosen


def _load_payload_bytes(raw: bytes, label: str) -> dict:
    with gzip.GzipFile(fileobj=io.BytesIO(raw), mode='rb') as g:
        x = safe_load_builtin_pickle(g)
    if x.get('kind') != 'BQG_MIXED_MASTER_KEY_METADATA' or int(x.get('schema_version', -1)) != 1:
        raise RuntimeError(f'{label}: wrong keymeta kind/schema')
    return x


def load_payloads_from_zip(zpath: Path) -> list[dict]:
    with zipfile.ZipFile(zpath) as z:
        names = []
        for info in z.infolist():
            normalized = info.filename.replace('\\\\', '/')
            parts = Path(normalized).parts
            mode = (info.external_attr >> 16) & 0xFFFF
            if (Path(normalized).is_absolute() or '..' in parts
                    or stat.S_ISLNK(mode) or stat.S_ISCHR(mode)
                    or stat.S_ISBLK(mode) or stat.S_ISFIFO(mode)):
                raise RuntimeError(f'{zpath}: unsafe ZIP member {info.filename!r}')
            if info.filename.endswith('.pkl.gz'):
                names.append(info.filename)
        names.sort()
        if not names:
            raise RuntimeError(f'{zpath}: expected at least one .pkl.gz')
        return [_load_payload_bytes(z.read(name), f'{zpath}:{name}') for name in names]


def load_payload_from_zip(zpath: Path) -> dict:
    payloads = load_payloads_from_zip(zpath)
    if len(payloads) != 1:
        raise RuntimeError(f'{zpath}: expected one .pkl.gz, got {len(payloads)}')
    return payloads[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', required=True)
    ap.add_argument('--run-id', action='append', type=int, required=True)
    ap.add_argument('--token', required=True)
    ap.add_argument('--work', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)

    chosen = choose_latest_keymeta(a.repo, a.run_id, a.token)
    payloads = []
    provenance = {}
    for sid in sorted(chosen):
        art = chosen[sid]
        zp = a.work / f'shard-{sid}.zip'
        size, actual_digest = download(art['archive_download_url'], a.token, zp)
        expected_digest = str(art.get('digest') or '')
        if not expected_digest.startswith('sha256:') or actual_digest != expected_digest.split(':', 1)[1].lower():
            raise RuntimeError(f"artifact {art['name']} SHA-256 mismatch or missing digest")
        if size != int(art.get('size_in_bytes') or -1):
            raise RuntimeError(f"artifact {art['name']} size mismatch")
        x = load_payload_from_zip(zp)
        if int(x['shard']) != sid:
            raise RuntimeError(f'artifact shard/name mismatch {sid} vs {x["shard"]}')
        payloads.append(x)
        provenance[str(sid)] = {
            'keymeta_artifact_id': art['id'],
            'keymeta_artifact_name': art['name'],
            'keymeta_artifact_digest': art.get('digest'),
            'keymeta_created_at': art.get('created_at'),
            'keymeta_run_id': art['source_keymeta_run_id'],
            'source_raw_run_id': x.get('source_run_id'),
            'source_raw_artifact_id': x.get('source_artifact_id'),
            'source_raw_artifact_digest': x.get('source_artifact_digest'),
            'source_tar_sha256': x.get('source_tar_sha256'),
        }
        zp.unlink(missing_ok=True)

    coverage = S.validate_keymeta_coverage(
        payloads,
        irrep='32',
        expected_shards=112,
        target_blocks=2755,
        target_columns=130903,
    )
    result = {
        'schema_version': 1,
        'kind': 'BQG_DEPTH6_32_ACTUAL_Q_KEYMETA_INVENTORY',
        'source_keymeta_runs': sorted(set(a.run_id)),
        **coverage,
        'provenance_by_shard': provenance,
        'structural_recompute_performed': False,
        'schedule_started': False,
        'next_action': (
            'GLOBAL_SCHEDULE_MAY_START'
            if coverage['schedule_allowed']
            else 'RECOVER_MISSING_KEYMETA_SHARDS_ONLY'
        ),
        'claim_boundary': (
            'Inventory of persisted actual-q metadata only. Partial coverage is not a rank witness. '
            'Global schedule/SVD remains forbidden until schedule_allowed=true.'
        ),
    }
    a.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({
        k: result[k] for k in (
            'status','recovered_shards','missing_shards','recovered_blocks',
            'recovered_columns','target_blocks','target_columns','schedule_allowed',
            'next_action')
    }, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
