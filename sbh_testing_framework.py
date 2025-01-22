#!/usr/bin/env python3

"""
improved_framework.py - Enhanced Testing Framework for SBH Solvers with CLI Support
"""

import os
import sys
import sqlite3
import json
import time
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path


class SBHTestingFramework:
    def __init__(self, db_path, test_dir):
        self.db_path = db_path
        self.test_dir = test_dir
        self.solver_path = None  # Inicjalizacja ścieżki solvera
        self._init_database()

    def _init_database(self):
        """Initialize SQLite database with required schema"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS solver_versions (
                    version_id INTEGER PRIMARY KEY,
                    version_name TEXT UNIQUE,
                    description TEXT,
                    created_at DATETIME
                )
            """)

            conn.execute("""
                CREATE TABLE IF NOT EXISTS test_cases (
                    test_id INTEGER PRIMARY KEY,
                    difficulty TEXT,
                    instance_path TEXT,
                    created_at DATETIME
                )
            """)

            conn.execute("""
                CREATE TABLE IF NOT EXISTS test_results (
                    result_id INTEGER PRIMARY KEY,
                    version_id INTEGER,
                    test_id INTEGER,
                    levenshtein_distance INTEGER,
                    execution_time_ms INTEGER,
                    original_dna TEXT,
                    reconstructed_dna TEXT,
                    tested_at DATETIME,
                    FOREIGN KEY (version_id) REFERENCES solver_versions(version_id),
                    FOREIGN KEY (test_id) REFERENCES test_cases(test_id)
                )
            """)
            conn.commit()

    def generate_test_cases(self):
        """Generate test cases using the C++ program's generate_instance mode"""
        if self.test_dir.exists() and any(self.test_dir.iterdir()):
            print("Przypadki testowe już istnieją, pomijanie generowania.")
            return

        os.makedirs(self.test_dir, exist_ok=True)

        # Definiowanie konfiguracji testów dla różnych poziomów trudności
        test_configs = {
            "easy": {
                "counts": 20,  # Zamiast 5
                "params": [
                    {"n": 300, "k": 7, "dk": 0, "ln": 10, "lp": 10},
                    {"n": 325, "k": 7, "dk": 0, "ln": 11, "lp": 11},
                    {"n": 350, "k": 8, "dk": 0, "ln": 12, "lp": 12},
                    {"n": 375, "k": 8, "dk": 0, "ln": 13, "lp": 13}
                ]
            },
            "medium": {
                "counts": 15,  # Zamiast 5
                "params": [
                    {"n": 400, "k": 8, "dk": 1, "ln": 15, "lp": 15},
                    {"n": 425, "k": 8, "dk": 1, "ln": 16, "lp": 16},
                    {"n": 450, "k": 9, "dk": 1, "ln": 18, "lp": 18},
                    {"n": 475, "k": 9, "dk": 1, "ln": 19, "lp": 19}
                ]
            },
            "hard": {
                "counts": 10,  # Zamiast 5
                "params": [
                    {"n": 500, "k": 9, "dk": 2, "ln": 20, "lp": 20},
                    {"n": 550, "k": 9, "dk": 2, "ln": 22, "lp": 22},
                    {"n": 600, "k": 10, "dk": 2, "ln": 25, "lp": 25},
                    {"n": 650, "k": 10, "dk": 2, "ln": 27, "lp": 27}
                ]
            }
        }

        # Generowanie przypadków testowych za pomocą programu C++
        with sqlite3.connect(self.db_path) as conn:
            for difficulty, config in test_configs.items():
                diff_dir = self.test_dir / difficulty
                diff_dir.mkdir(parents=True, exist_ok=True)

                for i in range(config["counts"]):
                    for params in config["params"]:
                        test_subdir = diff_dir / f"test_{i + 1}"
                        test_subdir.mkdir(parents=True, exist_ok=True)

                        instance_path = test_subdir / "instance.txt"

                        # Generowanie instancji za pomocą programu C++
                        cmd = (
                            f"{self.solver_path} generate_instance "
                            f"-n {params['n']} -k {params['k']} -dk {params['dk']} "
                            f"-ln {params['ln']} -lp {params['lp']} "
                            f"-o {instance_path}"
                        )

                        try:
                            subprocess.run(cmd, shell=True, check=True)
                        except subprocess.CalledProcessError as e:
                            print(f"Błąd podczas generowania instancji {instance_path}: {e}")
                            continue

                        # Zapisanie informacji o przypadku testowym do bazy danych
                        conn.execute("""
                            INSERT INTO test_cases (difficulty, instance_path, created_at)
                            VALUES (?, ?, ?)
                        """, (difficulty, str(instance_path), datetime.now()))

            conn.commit()

    def run_single_test(self, args):
        """Run solver on a single test case"""
        solver_path, instance_path, version_id = args
        instance_path = Path(instance_path)
        results_path = instance_path.with_name("results.txt")

        start_time = time.time()

        # Uruchomienie solvera w trybie test_instance
        cmd = f"{solver_path} test_instance -i {instance_path} -o {results_path}"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(f"Błąd podczas uruchamiania solvera na {instance_path}")
            return None

        execution_time = int((time.time() - start_time) * 1000)  # w milisekundach

        # Parsowanie wyników
        try:
            with open(results_path, 'r') as f:
                lines = f.readlines()
                original_dna = lines[0].split(": ")[1].strip()
                reconstructed_dna = lines[1].split(": ")[1].strip()
                levenshtein = int(lines[2].split(": ")[1].strip())

            return {
                'instance_path': str(instance_path),
                'execution_time': execution_time,
                'original_dna': original_dna,
                'reconstructed_dna': reconstructed_dna,
                'levenshtein_distance': levenshtein,
                'version_id': version_id
            }
        except Exception as e:
            print(f"Błąd podczas przetwarzania wyników dla {instance_path}: {str(e)}")
            return None

    def test_solver_version(self, solver_path, version_name, description, num_processes=None):
        """Test solver version using multiple processes"""
        if num_processes is None:
            num_processes = os.cpu_count()

        self.solver_path = solver_path

        # Rejestracja wersji solvera
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO solver_versions (version_name, description, created_at)
                VALUES (?, ?, ?)
                ON CONFLICT(version_name) DO UPDATE SET
                    description=excluded.description,
                    created_at=excluded.created_at
            """, (version_name, description, datetime.now()))
            conn.commit()
            cursor.execute("SELECT version_id FROM solver_versions WHERE version_name = ?", (version_name,))
            version_id = cursor.fetchone()[0]

        # Zbieranie wszystkich przypadków testowych
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT instance_path FROM test_cases")
            test_cases = [(self.solver_path, row[0], version_id) for row in cursor.fetchall()]

        # Uruchamianie testów równolegle
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            results = list(executor.map(self.run_single_test, test_cases))

        # Zapisanie wyników do bazy danych
        with sqlite3.connect(self.db_path) as conn:
            for result in results:
                if result:
                    conn.execute("""
                        INSERT INTO test_results 
                        (version_id, test_id, levenshtein_distance, execution_time_ms, 
                         original_dna, reconstructed_dna, tested_at)
                        VALUES (
                            ?, 
                            (SELECT test_id FROM test_cases WHERE instance_path = ?),
                            ?, ?, ?, ?, ?
                        )
                    """, (
                        result['version_id'],
                        result['instance_path'],
                        result['levenshtein_distance'],
                        result['execution_time'],
                        result['original_dna'],
                        result['reconstructed_dna'],
                        datetime.now()
                    ))
            conn.commit()

    def display_rankings(self):
        """Display rankings of all solver versions"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT 
                    v.version_name,
                    v.description,
                    COUNT(r.result_id) as total_tests,
                    AVG(r.execution_time_ms) as avg_execution_time,
                    AVG(CASE WHEN r.levenshtein_distance = 0 THEN 1.0 ELSE 0.0 END) as perfect_solutions,
                    AVG(r.levenshtein_distance) as avg_levenshtein
                FROM solver_versions v
                LEFT JOIN test_results r ON v.version_id = r.version_id
                GROUP BY v.version_id
                ORDER BY perfect_solutions DESC, avg_levenshtein ASC, avg_execution_time ASC
            """)

            results = cursor.fetchall()

            print("\nSOLVER RANKINGS")
            print("=" * 100)
            print(f"{'Version':<30} {'Tests':<8} {'Avg Time':<12} {'Perfect':<10} {'Avg Levenshtein':<15}")
            print("-" * 100)

            for row in results:
                version, desc, tests, avg_time, perfect, avg_lev = row
                perfect_percentage = (perfect * 100) if perfect is not None else 0
                avg_time = avg_time if avg_time is not None else 0
                avg_lev = avg_lev if avg_lev is not None else float('inf')
                print(f"{version:<30} {tests:>7} {avg_time:>11.0f}ms {perfect_percentage:>9.1f}% {avg_lev:>14.2f}")


def main():
    parser = argparse.ArgumentParser(description="SBH Testing Framework")

    # Określenie katalogu, w którym znajduje się skrypt
    script_dir = Path(__file__).parent

    # Poprawiona domyślna ścieżka do solvera w katalogu build wewnątrz katalogu projektu
    build_dir = script_dir / 'build'
    default_solver = build_dir / 'optymalizacja_kombinatoryczna'

    # Domyślna ścieżka do bazy danych i test_cases w build/
    default_db = build_dir / 'results.db'
    default_test_dir = build_dir / 'test_cases'

    parser.add_argument('--generate', action='store_true', help='Generate test cases (if not exist)')
    parser.add_argument('--test', action='store_true', help='Run tests on solver')
    parser.add_argument('--solver', default=str(default_solver), help='Path to solver executable')
    parser.add_argument('--version', help='Solver version name')
    parser.add_argument('--desc', help='Solver version description')
    parser.add_argument('--processes', type=int, help='Number of parallel processes')
    parser.add_argument('--rankings', action='store_true', help='Display solver rankings')

    args = parser.parse_args()

    # Wyświetlenie domyślnej ścieżki do solvera dla debugowania
    print(f"Domyślna ścieżka do solvera: {args.solver}")

    # Inicjalizacja frameworka z poprawionymi ścieżkami
    framework = SBHTestingFramework(db_path=str(default_db), test_dir=Path(args.solver).parent / 'test_cases')

    # Ustawienie poprawnej ścieżki do test_cases w build/
    framework.test_dir = build_dir / 'test_cases'

    if args.generate:
        # Ustawienie ścieżki solvera przed generowaniem testów
        framework.solver_path = args.solver
        framework.generate_test_cases()

    if args.test:
        if not args.version:
            print("Error: --version required for testing")
            sys.exit(1)
        framework.test_solver_version(args.solver, args.version,
                                      args.desc or "", args.processes)

    if args.rankings:
        framework.display_rankings()


if __name__ == "__main__":
    main()
