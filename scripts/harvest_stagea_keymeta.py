#!/usr/bin/env python3
"""Sequentially salvage tiny key metadata from a Stage-A Actions run.

The script intentionally never downloads all raw artifacts at once. It lists
matching shard artifacts, downloads ONE ZIP, extracts its one raw tar, runs the
persisted-map key metadata extractor, then deletes the multi-GB temporary data
before moving to the next shard.

It accepts partial runs: every successfully harvested shard is preserved and a
manifest records missing/expired/failed shards. Full numerical scheduling still
requires exact 112/112 metadata coverage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'proof_tools'))
from safe_archive import safe_extract_zip


API_HEADERS = {
    'Accept': 'application/vnd.github+json',
    'X-GitHub-Api-Version': '2022-11-28',
    'User-Agent': 'bqg-stagea-keymeta-harvester',
}


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def api_json(url: str, token: str):
    req = urllib.request.Request(url, headers={
        **API_HEADERS,
        'Authorization': f'Bearer {token}',
    })
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def list_artifacts(repo: str, run_id: str, token: str):
    out = []
    page = 1
    while True:
        x = api_json(
            f'https://api.github.com/repos/{repo}/actions/runs/{run_id}/artifacts?per_page=100&page={page}',
            token,
        )
        batch = x.get('artifacts', [])
        out.extend(batch)
        if len(batch) < 100:
            break
        page += 1
    return out


def download(url: str, token: str, dst: Path):
    """Download a GitHub Actions artifact without leaking GitHub auth to blob storage.

    GitHub's artifact endpoint returns a signed redirect URL on a storage host.
    Sending the GitHub Authorization header to that host can invalidate the
    signed request (401). Resolve the API redirect with auth, then download the
    signed URL without Authorization.
    """
    req = urllib.request.Request(url, headers={
        **API_HEADERS,
        'Authorization': f'Bearer {token}',
    })
    opener = urllib.request.build_opener(_NoRedirect())
    try:
        r = opener.open(req, timeout=120)
    except urllib.error.HTTPError as exc:
        if exc.code not in (301, 302, 303, 307, 308):
            raise
        location = exc.headers.get('Location')
        if not location:
            raise RuntimeError(f'artifact redirect {exc.code} missing Location') from exc
        if urllib.parse.urlparse(location).scheme.lower() != 'https':
            raise RuntimeError('refusing non-HTTPS signed artifact URL')
        signed_req = urllib.request.Request(location, headers={
            'User-Agent': API_HEADERS['User-Agent'],
        })
        with urllib.request.urlopen(signed_req, timeout=180) as r2, dst.open('wb') as f:
            shutil.copyfileobj(r2, f, length=8 * 1024 * 1024)
        return
    else:
        with r, dst.open('wb') as f:
            shutil.copyfileobj(r, f, length=8 * 1024 * 1024)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', required=True)
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--token-env', default='GH_TOKEN')
    ap.add_argument('--prefix', default='bqg-mixed-32-shard-')
    ap.add_argument('--expected-shards', type=int, default=112)
    ap.add_argument('--extractor', type=Path, required=True)
    ap.add_argument('--out-dir', type=Path, required=True)
    ap.add_argument('--engine-git-blob-sha', required=True)
    a = ap.parse_args()
    token = os.environ.get(a.token_env)
    if not token:
        raise RuntimeError(f'missing token env {a.token_env}')
    a.out_dir.mkdir(parents=True, exist_ok=True)

    rx = re.compile(r'^' + re.escape(a.prefix) + r'(\d+)$')
    artifacts = list_artifacts(a.repo, a.run_id, token)
    selected = {}
    duplicate = {}
    for art in artifacts:
        m = rx.match(str(art.get('name', '')))
        if not m:
            continue
        shard = int(m.group(1))
        if shard in selected:
            duplicate.setdefault(shard, []).append(int(art['id']))
            old = selected[shard]
            if str(art.get('created_at', '')) > str(old.get('created_at', '')):
                selected[shard] = art
        else:
            selected[shard] = art

    manifest = {
        'schema_version': 2,
        'kind': 'BQG_STAGEA_KEYMETA_SALVAGE_MANIFEST',
        'repo': a.repo,
        'source_run_id': str(a.run_id),
        'prefix': a.prefix,
        'expected_shards': a.expected_shards,
        'discovered_shards': sorted(selected),
        'duplicates': duplicate,
        'harvested': {},
        'errors': {},
    }

    for shard in sorted(selected):
        art = selected[shard]
        if art.get('expired'):
            manifest['errors'][str(shard)] = 'artifact_expired_before_harvest'
            continue
        tmp_root = Path(tempfile.mkdtemp(prefix=f'bqg-keymeta-{shard}-'))
        try:
            zpath = tmp_root / f'shard-{shard}.zip'
            download(str(art['archive_download_url']), token, zpath)
            zdir = tmp_root / 'zip'
            zdir.mkdir()
            expected_digest = str(art.get('digest') or '')
            if not expected_digest.startswith('sha256:'):
                raise RuntimeError('artifact metadata has no SHA-256 digest')
            actual_digest = sha256_file(zpath)
            if actual_digest != expected_digest.split(':', 1)[1].lower():
                raise RuntimeError('artifact ZIP SHA-256 mismatch')
            if zpath.stat().st_size != int(art.get('size_in_bytes') or -1):
                raise RuntimeError('artifact ZIP size mismatch')
            with zipfile.ZipFile(zpath) as zf:
                safe_extract_zip(zf, zdir)
            raw = zdir / f'{a.prefix}{shard}.tar.gz'
            if not raw.is_file():
                matches = list(zdir.rglob(f'{a.prefix}{shard}.tar.gz'))
                if len(matches) != 1:
                    raise RuntimeError(f'expected one raw tar for shard {shard}, got {len(matches)}')
                raw = matches[0]
            out_pkl = a.out_dir / f'{a.prefix}{shard}-keymeta.pkl.gz'
            out_json = a.out_dir / f'{a.prefix}{shard}-keymeta-stats.json'
            work = tmp_root / 'work'
            cmd = [
                sys.executable, str(a.extractor),
                '--src', str(raw), '--work', str(work),
                '--out', str(out_pkl), '--json', str(out_json),
                '--source-run-id', str(a.run_id),
                '--source-artifact-id', str(art['id']),
                '--source-artifact-digest', str(art.get('digest') or ''),
                '--engine-git-blob-sha', a.engine_git_blob_sha,
            ]
            child_env = os.environ.copy()
            child_env.pop(a.token_env, None)
            child_env.pop('GITHUB_TOKEN', None)
            subprocess.run(cmd, check=True, env=child_env)
            stats = json.loads(out_json.read_text())
            manifest['harvested'][str(shard)] = {
                'artifact_id': int(art['id']),
                'artifact_digest': art.get('digest'),
                'artifact_size_bytes': int(art.get('size_in_bytes') or 0),
                'artifact_expires_at': art.get('expires_at'),
                'blocks': int(stats['stats']['blocks']),
                'columns': int(stats['stats']['columns']),
                'q_keys': int(stats['stats']['q_keys']),
                'keymeta_bytes': int(stats['stats']['output_bytes']),
                'source_tar_sha256': stats['source_tar_sha256'],
            }
            print('HARVESTED', shard, manifest['harvested'][str(shard)], flush=True)
        except Exception as exc:
            manifest['errors'][str(shard)] = repr(exc)
            print('HARVEST_ERROR', shard, repr(exc), flush=True)
        finally:
            shutil.rmtree(tmp_root, ignore_errors=True)

    expected = set(range(a.expected_shards))
    harvested = {int(x) for x in manifest['harvested']}
    discovered = set(selected)
    manifest['harvested_shards'] = sorted(harvested)
    manifest['missing_shards'] = sorted(expected - harvested)
    manifest['harvested_count'] = len(harvested)
    manifest['harvested_blocks'] = sum(x['blocks'] for x in manifest['harvested'].values())
    manifest['harvested_columns'] = sum(x['columns'] for x in manifest['harvested'].values())
    manifest['status'] = 'COMPLETE_112_OF_112' if harvested == expected else 'PARTIAL_REUSABLE_KEYMETA'
    manifest['claim_boundary'] = (
        'Salvage of actual persisted numerical master-map key metadata only. '
        'Partial coverage is reusable but cannot prove [3,2] numerical closure.'
    )
    manifest_path = a.out_dir / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: manifest[k] for k in (
        'status','harvested_count','harvested_blocks','harvested_columns','missing_shards')}, indent=2))

    # Partial salvage is allowed, but a run that discovers existing raw shards
    # and preserves none of them is an infrastructure failure, not success.
    if discovered and not harvested:
        print('SALVAGE_FATAL: discovered raw shards but harvested zero', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
