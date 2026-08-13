#!/usr/bin/env python3
"""Focused synthetic tests for radiative-table validation."""

import io
import unittest

from check_radcorr_tables import validate_table


HEADER = "Ebeam,nu,XSborn_unp,XSrad_unp\n"


class ValidateTableTests(unittest.TestCase):
    def test_rejects_nonuniform_finite_grid(self):
        content = HEADER + "\n".join(
            f"10380,{nu},1.0,1.0" for nu in (6000, 6001, 6003, 6004)
        ) + "\n"
        with self.assertRaisesRegex(ValueError, "nonuniform finite nu grid"):
            validate_table("synthetic.csv", io.StringIO(content))

    def test_accepts_one_mev_grid_ending_at_9910(self):
        content = HEADER + "\n".join(
            f"10380,{nu},1.0,2.0" for nu in (9908, 9909, 9910)
        ) + "\n"
        summary = validate_table("synthetic.csv", io.StringIO(content))
        self.assertEqual(summary, (3, 9908.0, 9910.0, 3, 1.0))


if __name__ == "__main__":
    unittest.main()
