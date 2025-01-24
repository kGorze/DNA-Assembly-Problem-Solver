#!/usr/bin/env python3

import os
import sys
import sqlite3
import time
import argparse
import subprocess
import multiprocessing
import threading
import logging
from logging.handlers import RotatingFileHandler
import signal
import re
import curses
from datetime import datetime
from pathlib import Path, PureWindowsPath
from typing import Optional, Dict, Any
from display_manager import UnifiedDisplayManager

# Global reference to the framework instance
framework_instance = None

def setup_logger(debug_mode=False):
    # Konfiguracja głównego loggera
    logger = logging.getLogger()
    
    # Ustaw poziom logowania w zależności od trybu
    if debug_mode:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)
    
    # Usuń wszystkie istniejące handlery
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Utworzenie nowego handlera pliku z trybem 'w' (nadpisywanie)
    file_handler = logging.FileHandler('sbh_framework.log', mode='w', encoding='utf-8')
    
    if debug_mode:
        file_handler.setLevel(logging.DEBUG)
    else:
        file_handler.setLevel(logging.WARNING)
    
    # Format logów
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Dodanie handlera do loggera
    logger.addHandler(file_handler)
    
    return logger

# Użyj na początku programu
logger = setup_logger()

def cleanup_and_exit(signum, frame):
    if framework_instance and framework_instance.display_manager:
        framework_instance.display_manager.cleanup()
    try:
        curses.endwin()
    except:
        pass
    sys.exit(0)

signal.signal(signal.SIGINT, cleanup_and_exit)

def signal_handler(signum: int, frame: Any) -> None:
    global framework_instance
    logging.info(f"Received signal {signum}. Initiating cleanup.")
    if framework_instance and framework_instance.display_manager:
        framework_instance.display_manager.cleanup()
    print("\nReceived interrupt signal. Cleaning up...")
    sys.exit(1)

signal.signal(signal.SIGINT, signal_handler)

# Nowa regex do wykrywania pol rozbudowanego PROGRESS_UPDATE:
progress_pattern_expanded = re.compile(
    r"^PROGRESS_UPDATE:(\d+):([0-9\.]+):(.+?):([0-9\.]+):([0-9\.]+):([0-9\.]+):([0-9\.]+)$"
)
# Tłumaczenie:
# 1 -> pid
# 2 -> progress
# 3 -> status
# 4 -> bestFitness
# 5 -> coverage
# 6 -> edgeScore
# 7 -> theoreticalMax

def to_windows_path(path: Path) -> str:
    """Convert a path to Windows format."""
    path_str = str(path)
    # Jeśli jesteśmy w WSL, używamy oryginalnej ścieżki
    if path_str.startswith('/mnt/'):
        return path_str
    if sys.platform == 'win32':
        return str(PureWindowsPath(path))
    return str(path)

def read_config(config_path: Path) -> Dict[str, str]:
    """Czyta parametry z pliku konfiguracyjnego."""
    config = {}
    try:
        logging.warning(f"Reading config from: {config_path}")
        if not config_path.exists():
            logging.error(f"Config file not found at: {config_path}")
            return config
            
        logging.warning(f"Config file exists, size: {config_path.stat().st_size} bytes")
        with open(config_path, 'r') as f:
            content = f.read()
            logging.warning(f"Raw config content:\n{content}")
            
            for line in content.splitlines():
                line = line.strip()
                if line and not line.startswith('#'):
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        if value:
                            config[key] = value
                            logging.warning(f"Config: parsed {key} = {value}")
                        else:
                            logging.error(f"Config: empty value for key {key}")
    except Exception as e:
        logging.error(f"Error reading config file: {e}")
        import traceback
        logging.error(traceback.format_exc())
    
    if not config:
        logging.error("No valid configuration parameters found!")
    else:
        logging.warning(f"Successfully parsed {len(config)} parameters")
    return config

def solver_subprocess(
        solver_path: str,
        instance_path: Path,
        version_id: int,
        test_id: int,
        process_id: int,
        db_path: str,
        queue: multiprocessing.Queue
) -> None:
    start_time = time.time()
    
    # Określ trudność na podstawie ścieżki
    difficulty = "Unknown"
    path_str = str(instance_path).lower()
    if "easy" in path_str:
        difficulty = "Easy"
    elif "medium" in path_str:
        difficulty = "Medium"
    elif "hard" in path_str:
        difficulty = "Hard"
    
    logging.warning(f"Process {process_id}: Starting solver for test_id {test_id}, difficulty: {difficulty}")

    try:
        # Konwertuj ścieżki
        solver_path = to_windows_path(Path(solver_path))
        instance_path_str = to_windows_path(instance_path)
        results_path_str = to_windows_path(instance_path.with_name("results.txt"))
        
        # Wczytaj konfigurację
        config_path = Path(solver_path).parent / "config.cfg"
        logging.warning(f"Process {process_id}: Looking for config at: {config_path} (exists: {config_path.exists()})")
        if not config_path.exists():
            logging.error(f"Process {process_id}: Config file not found!")
            alternative_path = Path.cwd() / "build" / "config.cfg"
            if alternative_path.exists():
                logging.warning(f"Process {process_id}: Found alternative config at: {alternative_path}")
                config_path = alternative_path
            
        config = read_config(config_path)
        if not config:
            logging.error(f"Process {process_id}: No configuration loaded!")
        
        # Przygotuj parametry z konfiguracji
        config_params = []
        instance_params = ['k', 'deltaK', 'lNeg', 'lPoz', 'repAllowed', 'probablePositive']
        for key, value in config.items():
            if key not in instance_params:  # parametry instancji
                # Użyj formatu -param value zamiast --param=value
                param = f"-{key}"
                config_params.extend([param, value])
                logging.warning(f"Process {process_id}: Adding parameter: {param} {value}")
            else:
                logging.warning(f"Process {process_id}: Skipping instance parameter: {key}={value}")

        cmd = [
            solver_path,
            "test_instance",
            "-i", instance_path_str,
            "-o", results_path_str,
            "-pid", str(process_id),
            *config_params  # dodaj parametry z konfiguracji
        ]

        logging.warning(f"Process {process_id}: Full command: {' '.join(cmd)}")
        
        # Sprawdź czy pliki istnieją przed uruchomieniem
        logging.warning(f"Process {process_id}: Checking files:")
        logging.warning(f"  Solver path: {solver_path}")
        logging.warning(f"  Input path: {instance_path_str}")
        logging.warning(f"  Output path: {results_path_str}")
        logging.warning(f"  Solver exists: {Path(solver_path).exists()}")
        logging.warning(f"  Input exists: {Path(instance_path_str).exists()}")
        logging.warning(f"  Output dir exists: {Path(results_path_str).parent.exists()}")
        
        # Uruchom solver z przekierowaniem stderr do stdout
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1,
            shell=False  # Nie używaj powłoki
        )

        # Wysyłamy informację o trudności od razu po starcie
        queue.put(("current_test", process_id, f"CURRENT_TEST:{process_id}:{difficulty}:{instance_path}"))
        queue.put(("started", process_id))

        # Zbieraj wszystkie linie wyjścia
        all_output = []
        for line in proc.stdout:
            line = line.strip()
            if line:
                all_output.append(line)
                logging.warning(f"Process {process_id} output: {line}")
            if line.startswith("PROGRESS_UPDATE:"):
                match = progress_pattern_expanded.match(line)
                if match:
                    try:
                        pid = int(match.group(1))
                        progress_val = float(match.group(2))
                        status = match.group(3)
                        best_fitness = float(match.group(4))
                        coverage_val = float(match.group(5))
                        edge_score_val = float(match.group(6))
                        theoretical_max = float(match.group(7))

                        queue.put(("progress_expanded",
                                   pid,
                                   progress_val,
                                   status,
                                   best_fitness,
                                   coverage_val,
                                   edge_score_val,
                                   theoretical_max))
                    except (ValueError, IndexError) as e:
                        error_msg = f"Process {process_id}: Error converting expanded progress update: {e}"
                        queue.put(("output", process_id, error_msg, True))
                        logging.error(error_msg)
                else:
                    error_msg = f"Process {process_id}: Unrecognized progress format: {line}"
                    queue.put(("output", process_id, error_msg, True))
                    logging.error(error_msg)
            else:
                queue.put(("output", process_id, line, False))

        proc.wait()
        exit_code = proc.returncode
        logging.info(f"Process {process_id}: Solver exited with code {exit_code}.")

        if exit_code != 0:
            error_msg = f"Process {process_id}: Solver exited with code {exit_code}"
            queue.put(("output", process_id, error_msg, True))
            logging.error(error_msg)

        # Wczytanie results.txt
        results_path = instance_path.with_name("results.txt")
        dist = -1
        original_dna = ""
        reconstructed_dna = ""

        try:
            if results_path.exists():
                with open(results_path, "r") as f:
                    lines = f.readlines()
                    if len(lines) >= 3:
                        original_dna = lines[0].split(": ", 1)[1].strip()
                        reconstructed_dna = lines[1].split(": ", 1)[1].strip()
                        dist = int(lines[2].split(": ", 1)[1].strip())
                        logging.info(f"Process {process_id}: Results - dist={dist}")
            else:
                logging.warning(f"Process {process_id}: Results file not found: {results_path}")
        except Exception as e:
            error_msg = f"Process {process_id}: Error reading results file: {e}"
            queue.put(("output", process_id, error_msg, True))
            logging.error(error_msg)

        exec_time_ms = int((time.time() - start_time) * 1000)
        logging.info(f"Process {process_id}: Execution time: {exec_time_ms} ms.")

        queue.put((
            "done",
            process_id,
            test_id,
            dist,
            original_dna,
            reconstructed_dna,
            exec_time_ms
        ))
        logging.info(f"Process {process_id}: Sent 'done' message to main process.")

    except Exception as e:
        error_msg = f"Process {process_id}: Error during solver run: {e}"
        queue.put(("output", process_id, error_msg, True))
        logging.exception(error_msg)

class SBHTestingFramework:
    def __init__(self, db_path: str, test_dir: str, use_ui: bool = True) -> None:
        self.db_path = db_path
        self.test_dir = Path(test_dir)
        self.solver_path: Optional[str] = None
        self.display_manager: Optional[UnifiedDisplayManager] = None
        self.use_ui = use_ui and sys.stdout.isatty()
        logging.info("Initializing SBH Testing Framework.")
        try:
            self._init_database()
            logging.info("Database initialized successfully.")
        except sqlite3.Error as e:
            logging.error(f"Database initialization failed: {e}")
            sys.exit(1)

    def _init_database(self) -> None:
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

    def generate_test_cases(self) -> None:
        if self.test_dir.exists() and any(self.test_dir.iterdir()):
            print("Test cases already exist, skipping generation.")
            logging.info("Test cases already exist, skipping generation.")
            return

        os.makedirs(self.test_dir, exist_ok=True)
        logging.info(f"Created test directory at {self.test_dir}.")

        test_configs: Dict[str, Dict[str, Any]] = {
            "easy": {
                "counts": 5,  # 5 testów dla każdego zestawu parametrów
                "params": [
                    {"n": 300, "k": 7, "dk": 0, "ln": 10, "lp": 10},
                    {"n": 325, "k": 7, "dk": 0, "ln": 11, "lp": 11},
                    {"n": 350, "k": 8, "dk": 0, "ln": 12, "lp": 12},
                    {"n": 375, "k": 8, "dk": 0, "ln": 13, "lp": 13}
                ]
            },
            "medium": {
                "counts": 5,
                "params": [
                    {"n": 400, "k": 8, "dk": 1, "ln": 15, "lp": 15},
                    {"n": 425, "k": 8, "dk": 1, "ln": 16, "lp": 16},
                    {"n": 450, "k": 9, "dk": 1, "ln": 18, "lp": 18},
                    {"n": 475, "k": 9, "dk": 1, "ln": 19, "lp": 19}
                ]
            },
            "hard": {
                "counts": 5,
                "params": [
                    {"n": 500, "k": 9, "dk": 2, "ln": 20, "lp": 20},
                    {"n": 550, "k": 9, "dk": 2, "ln": 22, "lp": 22},
                    {"n": 600, "k": 10, "dk": 2, "ln": 25, "lp": 25},
                    {"n": 650, "k": 10, "dk": 2, "ln": 27, "lp": 27}
                ]
            }
        }

        total_test_count = 0
        with sqlite3.connect(self.db_path) as conn:
            for difficulty, config in test_configs.items():
                diff_dir = self.test_dir / difficulty
                diff_dir.mkdir(parents=True, exist_ok=True)
                logging.info(f"Creating directory for difficulty '{difficulty}' at {diff_dir}.")

                # Dla każdego zestawu parametrów
                for param_idx, params in enumerate(config["params"]):
                    # Generujemy counts testów z tymi parametrami
                    for i in range(config["counts"]):
                        test_subdir = diff_dir / f"test_{total_test_count + 1}"
                        test_subdir.mkdir(parents=True, exist_ok=True)

                        instance_path = test_subdir / "instance.txt"
                        cmd = (
                            f"{self.solver_path} generate_instance "
                            f"-n {params['n']} -k {params['k']} -dk {params['dk']} "
                            f"-ln {params['ln']} -lp {params['lp']} "
                            f"-o {instance_path}"
                        )
                        try:
                            logging.debug(f"Executing command: {cmd}")
                            subprocess.run(cmd, shell=True, check=True)
                            total_test_count += 1
                            logging.info(f"Generated instance at {instance_path} ({total_test_count} total)")
                        except subprocess.CalledProcessError as e:
                            error_msg = f"Error generating instance {instance_path}: {e}"
                            print(error_msg)
                            logging.error(error_msg)
                            continue

                        try:
                            conn.execute("""
                                INSERT INTO test_cases (difficulty, instance_path, created_at)
                                VALUES (?, ?, ?)
                            """, (difficulty, str(instance_path), datetime.now()))
                        except sqlite3.Error as e:
                            error_msg = f"Database error while inserting test case {instance_path}: {e}"
                            print(error_msg)
                            logging.error(error_msg)

            conn.commit()
            logging.info(f"All test cases generated and saved to the database. Total tests: {total_test_count}")
            print(f"Generated {total_test_count} test cases in total.")

    def test_solver_version(
            self,
            solver_path: str,
            version_name: str,
            description: str,
            num_processes: Optional[int] = None
    ) -> None:
        logging.info(f"Starting tests for solver version: {version_name}")
        try:
            self.solver_path = solver_path
            logging.info(f"Registering solver version '{version_name}' with description '{description}'.")
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute("""
                    INSERT OR IGNORE INTO solver_versions (version_name, description, created_at)
                    VALUES (?, ?, ?)
                """, (version_name, description, datetime.now()))
                conn.commit()

                cursor.execute("SELECT version_id FROM solver_versions WHERE version_name = ?", (version_name,))
                row = cursor.fetchone()
                if not row:
                    message = "Error: version not found in DB after insertion."
                    if self.use_ui and self.display_manager:
                        self.display_manager.add_output(message, error=True)
                    else:
                        logging.error(message)
                    return
                version_id = row[0]
                logging.info(f"Solver version '{version_name}' has version_id {version_id}.")

            with sqlite3.connect(self.db_path) as conn:
                c = conn.cursor()
                # Pobierz tylko najnowsze testy dla każdej trudności
                c.execute("""
                    WITH RankedTests AS (
                        SELECT 
                            test_id,
                            instance_path,
                            difficulty,
                            created_at,
                            ROW_NUMBER() OVER (PARTITION BY difficulty ORDER BY created_at DESC) as rn
                        FROM test_cases
                    )
                    SELECT test_id, instance_path, difficulty
                    FROM RankedTests
                    WHERE rn <= 20  -- 20 testów dla każdej trudności
                    ORDER BY difficulty, created_at
                """)
                test_cases = c.fetchall()
                
                # Logowanie informacji o testach
                difficulties = {}
                for _, _, diff in test_cases:
                    difficulties[diff] = difficulties.get(diff, 0) + 1
                logging.warning(f"Found {len(test_cases)} test cases to run:")
                for diff, count in difficulties.items():
                    logging.warning(f"- {diff}: {count} tests")

            if not test_cases:
                message = "No test cases found. Please generate first."
                if self.use_ui and self.display_manager:
                    self.display_manager.add_output(message, error=True)
                else:
                    logging.error(message)
                return

            total_tests = len(test_cases)
            logging.warning(f"Total test cases to run: {total_tests}")

            cpu_count = os.cpu_count() or 1
            if num_processes is None:
                num_processes = cpu_count
            else:
                num_processes = min(num_processes, cpu_count, total_tests)
            logging.info(f"Using {num_processes} parallel processes out of {cpu_count} available CPU cores.")

            self.display_manager = UnifiedDisplayManager(
                total_tests=total_tests,
                num_processes=num_processes,
                use_ui=self.use_ui
            )
            logging.info("Initialized UnifiedDisplayManager.")

            original_add_output = self.display_manager.add_output
            def new_add_output(message: str, error: bool = False) -> None:
                if error:
                    logging.error(message)
                else:
                    logging.info(message)
                original_add_output(message, error)

            self.display_manager.add_output = new_add_output
            logging.info("Enhanced display_manager.add_output to include logging.")

            self.display_manager.add_output(
                f"Using {num_processes} parallel processes out of {cpu_count} available CPU cores.",
                error=False
            )

            manager = multiprocessing.Manager()
            queue = manager.Queue()
            logging.info("Created multiprocessing.Queue for inter-process communication.")

            tests_to_run = list(test_cases)
            running_processes: list[multiprocessing.Process] = []
            test_index = 0
            logging.info("Prepared list of tests to run.")

            for pid in range(num_processes):
                if test_index >= total_tests:
                    break
                test_id, inst_path, difficulty = tests_to_run[test_index]
                p = multiprocessing.Process(
                    target=solver_subprocess,
                    args=(
                        solver_path,
                        Path(inst_path),
                        version_id,
                        test_id,
                        pid,
                        self.db_path,
                        queue
                    )
                )
                p.start()
                running_processes.append(p)
                logging.info(f"Started process {pid} for test_id {test_id}, difficulty: {difficulty}")
                test_index += 1

            stop_event = threading.Event()

            def queue_listener() -> None:
                nonlocal running_processes, test_index, total_tests, tests_to_run
                finished_tests = 0

                while not stop_event.is_set() and finished_tests < total_tests:
                    try:
                        msg = queue.get(timeout=1)
                    except:
                        continue
                    if not msg:
                        continue

                    mtype = msg[0]

                    if mtype == "progress_expanded":
                        # ("progress_expanded", pid, progress_val, status, best_fitness, coverage_val, edge_score_val, theoretical_max)
                        _, pid, progress_val, stat, best_fit, cov_val, edge_val, tmax = msg
                        self.display_manager.update_process(
                            pid, progress_val, stat, best_fit,
                            coverage=cov_val, edge_score=edge_val, theoretical_max=tmax
                        )
                        logging.debug(f"Main process: Received expanded progress from PID {pid}: {progress_val}% - {stat}")

                    elif mtype == "output":
                        _, pid, line, is_error = msg
                        self.display_manager.add_output(line, error=is_error)

                    elif mtype == "current_test":
                        # ("current_test", process_id, line)
                        _, pid, line = msg
                        try:
                            # Parse CURRENT_TEST:procid:difficulty:filepath
                            parts = line.split(":")
                            if len(parts) >= 4:
                                filepath = parts[3]
                                # Extract difficulty from filepath
                                difficulty = 'Unknown'  # default value
                                if '/easy/' in filepath:
                                    difficulty = 'Easy'
                                elif '/medium/' in filepath:
                                    difficulty = 'Medium'
                                elif '/hard/' in filepath:
                                    difficulty = 'Hard'
                    
                                # Safely update process status
                                if self.display_manager:
                                    if pid in self.display_manager.process_status:
                                        self.display_manager.process_status[pid]["difficulty"] = difficulty
                                        self.display_manager.process_status[pid]["current_test_name"] = filepath
                                        # Add logging to debug
                                        logging.debug(f"Updated process {pid} with difficulty: {difficulty}")
                    
                            # Always log the line regardless of parsing success
                            if self.display_manager:
                                self.display_manager.add_output(line, error=False)
                    
                        except Exception as e:
                            error_msg = f"Error processing current_test message for pid {pid}: {str(e)}"
                            logging.error(error_msg)
                            if self.display_manager:
                                self.display_manager.add_output(error_msg, error=True)

                    elif mtype == "done":
                        _, pid, done_test_id, dist, orig, recon, exec_ms = msg
                        logging.info(f"Main process: 'done' from PID {pid} test_id {done_test_id}, dist={dist}, time={exec_ms}ms")

                        try:
                            with sqlite3.connect(self.db_path) as conn:
                                c = conn.cursor()
                                c.execute("""
                                    INSERT INTO test_results
                                    (version_id, test_id, levenshtein_distance,
                                     execution_time_ms, original_dna, reconstructed_dna, tested_at)
                                    VALUES (?, ?, ?, ?, ?, ?, ?)
                                """, (
                                    version_id,
                                    done_test_id,
                                    dist,
                                    exec_ms,
                                    orig,
                                    recon,
                                    datetime.now()
                                ))
                                conn.commit()
                            logging.info(f"Saved results for test_id {done_test_id} to database.")
                        except sqlite3.Error as e:
                            error_msg = f"Database error while saving results for test_id {done_test_id}: {e}"
                            self.display_manager.add_output(error_msg, error=True)
                            logging.error(error_msg)

                        self.display_manager.add_output(
                            f"Test {done_test_id} (process {pid}) finished. dist={dist}",
                            error=False
                        )
                        self.display_manager.increment_completed()
                        finished_tests += 1

                        if test_index < total_tests:
                            new_test_id, new_inst_path, new_difficulty = tests_to_run[test_index]
                            p = multiprocessing.Process(
                                target=solver_subprocess,
                                args=(
                                    solver_path,
                                    Path(new_inst_path),
                                    version_id,
                                    new_test_id,
                                    pid,
                                    self.db_path,
                                    queue
                                )
                            )
                            p.start()
                            running_processes.append(p)
                            logging.info(f"Started new process {pid} for test_id {new_test_id}, difficulty: {new_difficulty}")
                            test_index += 1

                        active_procs = [p for p in running_processes if p.is_alive()]

                    elif mtype == "started":
                        # "started", process_id
                        pass

                logging.info("Queue listener has finished processing all tests.")

            listener_thread = threading.Thread(target=queue_listener, daemon=True)
            listener_thread.start()
            logging.info("Started queue listener thread.")

            try:
                while listener_thread.is_alive():
                    time.sleep(0.1)
            except KeyboardInterrupt:
                if self.display_manager:
                    self.display_manager.add_output("Interrupted by user.", error=True)
                logging.warning("Main process: Interrupted by user via KeyboardInterrupt.")
                raise

        except Exception as e:
            logging.error(f"Error in test_solver_version: {e}")
            raise
        finally:
            stop_event.set()
            if 'listener_thread' in locals():
                listener_thread.join()
                logging.info("Main process: Listener thread has been stopped.")

            for p in running_processes:
                if p.is_alive():
                    p.terminate()
                    logging.info(f"Main process: Terminated process {p.pid}.")
            for p in running_processes:
                p.join()
                logging.info(f"Main process: Joined process {p.pid}.")

            if self.display_manager:
                self.display_manager.cleanup()
                logging.info("Main process: Cleaned up display manager.")
            logging.info("Main process: Exiting test_solver_version method.")

    def display_rankings(self) -> None:
        logging.info("Displaying solver rankings.")
        try:
            with sqlite3.connect(self.db_path) as conn:
                c = conn.cursor()
                c.execute("""
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
                results = c.fetchall()

                print("\nSOLVER RANKINGS")
                print("=" * 100)
                print(f"{'Version':<30} {'Tests':<8} {'Avg Time':<12} {'Perfect':<10} {'Avg Lev':<15}")
                print("-" * 100)

                for row in results:
                    version, desc, tests, avg_time, perfect, avg_lev = row
                    perfect_percentage = (perfect * 100) if perfect is not None else 0
                    avg_time = avg_time if avg_time is not None else 0
                    avg_lev = avg_lev if avg_lev is not None else float('inf')
                    print(f"{version:<30} {tests:>7} {avg_time:>11.0f}ms {perfect_percentage:>9.1f}% {avg_lev:>14.2f}")

                logging.info("Displayed solver rankings successfully.")
        except sqlite3.Error as e:
            error_msg = f"SBHTestingFramework: Database error while displaying rankings: {e}"
            print(error_msg)
            logging.error(error_msg)
        except Exception as e:
            error_msg = f"SBHTestingFramework: Error while displaying rankings: {e}"
            print(error_msg)
            logging.error(error_msg)

def main() -> None:
    parser = argparse.ArgumentParser(description="SBH Testing Framework (multiprocessing + curses)")

    script_dir = Path(__file__).parent
    build_dir = script_dir / 'build'
    default_solver = build_dir / 'optymalizacja_kombinatoryczna'
    default_db = build_dir / 'results.db'
    default_test_dir = build_dir / 'test_cases'

    parser.add_argument('--generate', action='store_true', help='Generate test cases (if not exist)')
    parser.add_argument('--test', action='store_true', help='Run tests on solver in parallel')
    parser.add_argument('--solver', default=str(default_solver), help='Path to solver executable')
    parser.add_argument('--version', help='Solver version name')
    parser.add_argument('--desc', help='Solver version description')
    parser.add_argument('--processes', type=int, help='Number of parallel processes')
    parser.add_argument('--rankings', action='store_true', help='Display solver rankings')
    parser.add_argument('--log', action='store_true', help='Use logging instead of UI')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')

    args = parser.parse_args()

    try:
        use_ui = not args.log and sys.stdout.isatty()
        # Inicjalizacja loggera z odpowiednim trybem
        logger = setup_logger(debug_mode=args.debug)
        logging.info(f"Determined use_ui: {use_ui}")

        framework = SBHTestingFramework(
            db_path=str(default_db),
            test_dir=str(default_test_dir),
            use_ui=use_ui
        )
        global framework_instance
        framework_instance = framework

        framework.solver_path = args.solver
        logging.info(f"Using solver at path: {args.solver}")

        if not os.path.isfile(framework.solver_path):
            error_msg = f"Solver executable not found at path: {framework.solver_path}"
            print(error_msg)
            logging.error(error_msg)
            sys.exit(1)
        logging.info(f"Validated solver path: {framework.solver_path}")

        if args.generate:
            logging.info("Generating test cases as per --generate flag.")
            framework.generate_test_cases()
            logging.info("Test cases generation completed.")

        if args.test:
            logging.info("Running tests as per --test flag.")
            if not args.version:
                error_msg = "Error: --version required for testing"
                print(error_msg)
                logging.error(error_msg)
                sys.exit(1)
            framework.test_solver_version(
                solver_path=args.solver,
                version_name=args.version,
                description=args.desc or "",
                num_processes=args.processes
            )
            logging.info("Testing completed.")

        if args.rankings:
            logging.info("Displaying solver rankings as per --rankings flag.")
            framework.display_rankings()

    except Exception as e:
        error_msg = f"Main: Critical error in main: {e}"
        print(error_msg)
        logging.exception(error_msg)
        sys.exit(1)


if __name__ == "__main__":
    main()
