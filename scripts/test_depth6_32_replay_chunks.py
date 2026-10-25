#!/usr/bin/env python3
import unittest

import depth6_32_replay_chunks as C


class ReplayChunkPlanningTests(unittest.TestCase):
    def setUp(self):
        self.records = [
            {'orbit_id': 1, 'proxy': 100, 'm': 1, 'coord_dim': 100},
            {'orbit_id': 2, 'proxy': 90, 'm': 2, 'coord_dim': 45},
            {'orbit_id': 3, 'proxy': 80, 'm': 4, 'coord_dim': 20},
            {'orbit_id': 4, 'proxy': 70, 'm': 5, 'coord_dim': 14},
            {'orbit_id': 5, 'proxy': 400, 'm': 10, 'coord_dim': 40},
        ]

    def test_excludes_completed_and_covers_every_missing_orbit_once(self):
        chunks = C.plan_missing_chunks(self.records, {2, 4}, target_proxy=180)
        ids = [oid for ch in chunks for oid in ch['orbit_ids']]
        self.assertEqual(sorted(ids), [1, 3, 5])
        self.assertEqual(len(ids), len(set(ids)))

    def test_plan_is_deterministic_under_input_reordering(self):
        a = C.plan_missing_chunks(self.records, {2}, target_proxy=180)
        b = C.plan_missing_chunks(list(reversed(self.records)), {2}, target_proxy=180)
        self.assertEqual(a, b)

    def test_chunk_proxy_is_capped_except_oversize_singleton(self):
        chunks = C.plan_missing_chunks(self.records, set(), target_proxy=180)
        for ch in chunks:
            if ch['proxy'] > 180:
                self.assertEqual(len(ch['orbit_ids']), 1)
                self.assertEqual(ch['orbit_ids'], [5])
            else:
                self.assertLessEqual(ch['proxy'], 180)

    def test_max_blocks_bounds_work_per_chunk(self):
        chunks = C.plan_missing_chunks(self.records, set(), target_proxy=180, max_blocks=2)
        ids = [oid for ch in chunks for oid in ch['orbit_ids']]
        self.assertEqual(sorted(ids), [1, 2, 3, 4, 5])
        self.assertEqual(len(ids), len(set(ids)))
        for ch in chunks:
            self.assertLessEqual(len(ch['orbit_ids']), 2)
            self.assertEqual(ch['blocks'], len(ch['orbit_ids']))

    def test_unknown_completed_orbit_fails_closed(self):
        with self.assertRaises(RuntimeError):
            C.plan_missing_chunks(self.records, {999}, target_proxy=180)

    def test_duplicate_orbit_fails_closed(self):
        bad = self.records + [dict(self.records[0])]
        with self.assertRaises(RuntimeError):
            C.plan_missing_chunks(bad, set(), target_proxy=180)


if __name__ == '__main__':
    unittest.main()
