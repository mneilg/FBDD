#!/usr/bin/env python3

"""
reporter.py – Combined CSV Result Reporter and Progress Tracker for VinaScreen
Handles adding results to the summary CSV, prints a real-time progress count,
and writes error details to a separate log file.
"""

import csv
import os
from datetime import datetime

class CsvProgressReporter:
    def __init__(self, csv_file, total_jobs, log_file='errors.log'):
        self.csv_file = csv_file
        self.log_file = log_file
        self.total_jobs = total_jobs
        self.done = 0
        self.kept = 0
        self.skipped = 0
        self._write_header()
        # Clear the log file at the start of each run
        open(self.log_file, 'w').close()

    def _write_header(self):
        write_header = not os.path.exists(self.csv_file)
        with open(self.csv_file, "a", newline="") as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow([
                    "receptor", "ligand", "affinity", "rmsd_lb", "status", "output_file", "timestamp"
                ])

    def add_result(self, receptor, ligand, affinity, rmsd_lb, success, output_file=None, error_msg=None):
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.done += 1
        if success:
            self.kept += 1
            status = "OK"
        else:
            self.skipped += 1
            status = "FAILED"
            # Write error details to the log file
            with open(self.log_file, "a") as logf:
                logf.write(f"[{now}] {receptor}+{ligand}: {error_msg}\n")
        with open(self.csv_file, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                receptor, ligand,
                affinity if affinity is not None else "",
                rmsd_lb if rmsd_lb is not None else "",
                status, output_file if output_file else "", now
            ])
        self.report_progress()

    def report_progress(self):
        print(f"Processed: {self.done}/{self.total_jobs} (kept: {self.kept}, skipped: {self.skipped})",
              end='\r' if self.done != self.total_jobs else '\n', flush=True)

# Example usage for testing
if __name__ == "__main__":
    reporter = CsvProgressReporter("vinascreen_results.csv", total_jobs=3)
    reporter.add_result("rec1", "lig1", -8.1, 0.0, True, "output/rec1_lig1.pdbqt")
    reporter.add_result("rec1", "lig2", None, None, False, None, "Vina failed for lig2")
    reporter.add_result("rec2", "lig1", -7.2, 0.2, True, "output/rec2_lig1.pdbqt")

