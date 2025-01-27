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
from datetime import datetime, timedelta
from pathlib import Path, PureWindowsPath
from typing import Optional, Dict, Any, List
from display_manager import UnifiedDisplayManager
import shutil
from contextlib import contextmanager
import queue

# Constants for timeouts and retries
PROCESS_TIMEOUT = 300  # 5 minutes
QUEUE_TIMEOUT = 1.0   # 1 second
DB_RETRY_COUNT = 3
GRACEFUL_SHUTDOWN_TIMEOUT = 5.0

class ProcessTimeoutError(Exception):
    pass

class DatabaseError(Exception):
    pass

@contextmanager
def database_connection(db_path: str, timeout: float = 20.0):
    """Context manager for database connections with retry logic"""
    for attempt in range(DB_RETRY_COUNT):
        try:
            conn = sqlite3.connect(db_path, timeout=timeout)
            conn.execute("PRAGMA journal_mode=WAL")  # Enable Write-Ahead Logging
            yield conn
            conn.commit()
            break
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e) and attempt < DB_RETRY_COUNT - 1:
                time.sleep(0.1 * (attempt + 1))  # Exponential backoff
                continue
            raise DatabaseError(f"Database error after {attempt + 1} attempts: {e}")
        finally:
            if 'conn' in locals():
                conn.close()

class ProcessManager:
    def __init__(self):
        self.processes: List[multiprocessing.Process] = []
        self._lock = threading.Lock()
        
    def add_process(self, process: multiprocessing.Process):
        with self._lock:
            self.processes.append(process)
            
    def remove_process(self, process: multiprocessing.Process):
        with self._lock:
            if process in self.processes:
                self.processes.remove(process)
                
    def terminate_all(self, timeout: float = GRACEFUL_SHUTDOWN_TIMEOUT):
        """Terminate all processes with timeout"""
        with self._lock:
            for p in self.processes:
                if p.is_alive():
                    p.terminate()
                    
            # Wait for graceful termination
            deadline = time.time() + timeout
            while time.time() < deadline and any(p.is_alive() for p in self.processes):
                time.sleep(0.1)
                
            # Force kill any remaining processes
            for p in self.processes:
                if p.is_alive():
                    try:
                        p.kill()  # Force kill
                    except Exception as e:
                        logging.error(f"Failed to kill process {p.pid}: {e}")
            
            # Final cleanup
            for p in self.processes:
                try:
                    p.join(timeout=0.5)
                except Exception as e:
                    logging.error(f"Failed to join process {p.pid}: {e}")
            
            self.processes.clear()

class SolverProcess:
    def __init__(self, solver_path: str, process_id: int, timeout: float = PROCESS_TIMEOUT):
        self.solver_path = solver_path
        self.process_id = process_id
        self.timeout = timeout
        self.process: Optional[subprocess.Popen] = None
        self.start_time: Optional[float] = None
        
    def start(self, cmd: List[str], cwd: Optional[str] = None) -> None:
        self.process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            bufsize=1,
            cwd=cwd
        )
        self.start_time = time.time()
        
    def check_timeout(self) -> bool:
        if self.start_time and time.time() - self.start_time > self.timeout:
            return True
        return False
        
    def terminate(self) -> None:
        if self.process:
            try:
                self.process.terminate()
                try:
                    self.process.wait(timeout=1.0)
                except subprocess.TimeoutExpired:
                    self.process.kill()
            except Exception as e:
                logging.error(f"Error terminating process {self.process_id}: {e}")

class QueueManager:
    def __init__(self, queue_size: int = 1000):
        self.queue = multiprocessing.Queue(maxsize=queue_size)
        self.stop_event = multiprocessing.Event()
        self._lock = threading.Lock()
        
    def put(self, item: Any, timeout: float = QUEUE_TIMEOUT) -> bool:
        try:
            self.queue.put(item, timeout=timeout)
            return True
        except queue.Full:
            logging.warning("Queue is full, message dropped")
            return False
            
    def get(self, timeout: float = QUEUE_TIMEOUT) -> Optional[Any]:
        try:
            return self.queue.get(timeout=timeout)
        except queue.Empty:
            return None
            
    def stop(self) -> None:
        self.stop_event.set()

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
    try:
        logger = setup_logger(True)  # Enable debug mode for subprocess
        logger.info(f"Process {process_id}: Starting solver subprocess")
        logger.info(f"Process {process_id}: Solver path: {solver_path}")
        logger.info(f"Process {process_id}: Instance path: {instance_path}")
        
        # Validate paths
        if not os.path.exists(solver_path):
            logger.error(f"Process {process_id}: Solver executable not found at {solver_path}")
            raise FileNotFoundError(f"Solver not found: {solver_path}")
            
        if not os.path.exists(instance_path):
            logger.error(f"Process {process_id}: Instance file not found at {instance_path}")
            raise FileNotFoundError(f"Instance not found: {instance_path}")
            
        # Create output directory if it doesn't exist
        output_dir = Path("output")
        output_dir.mkdir(exist_ok=True)
        
        # Prepare output file path
        output_file = output_dir / f"result_{process_id}.txt"
        logger.info(f"Process {process_id}: Output will be written to {output_file}")
        
        # Prepare config file path and check if it exists
        config_file = Path("config.cfg")
        if not config_file.exists():
            logger.error(f"Process {process_id}: Config file not found at {config_file}")
            raise FileNotFoundError(f"Config file not found: {config_file}")
            
        logger.info(f"Process {process_id}: Using config file: {config_file}")
        
        # Prepare command
        cmd = [
            solver_path,
            "test_instance",
            "-i", str(instance_path),
            "-o", str(output_file),
            "-cfg", str(config_file),
            "-pid", str(process_id)
        ]
        
        logger.info(f"Process {process_id}: Executing command: {' '.join(cmd)}")
        
        # Execute solver
        try:
            start_time = time.time()
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            
            # Read output in real-time
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    logger.debug(f"Process {process_id} output: {output.strip()}")
                    
            # Get return code and remaining output
            return_code = process.wait()
            _, stderr = process.communicate()
            
            execution_time = time.time() - start_time
            
            if return_code != 0:
                logger.error(f"Process {process_id}: Solver failed with return code {return_code}")
                logger.error(f"Process {process_id}: Error output: {stderr}")
                raise RuntimeError(f"Solver failed with return code {return_code}")
                
            logger.info(f"Process {process_id}: Solver completed successfully in {execution_time:.2f} seconds")
            
            # Check if output file was created
            if not output_file.exists():
                logger.error(f"Process {process_id}: Output file was not created")
                raise FileNotFoundError(f"Output file not created: {output_file}")
                
            # Read results
            with open(output_file, 'r') as f:
                result = f.read()
                logger.info(f"Process {process_id}: Results: {result}")
                
            # Send results through queue
            queue.put({
                'process_id': process_id,
                'version_id': version_id,
                'test_id': test_id,
                'status': 'success',
                'execution_time': execution_time,
                'result': result
            })
            
        except subprocess.SubprocessError as e:
            logger.error(f"Process {process_id}: Failed to execute solver: {str(e)}")
            raise
            
    except Exception as e:
        logger.error(f"Process {process_id}: Fatal error: {str(e)}")
        queue.put({
            'process_id': process_id,
            'version_id': version_id,
            'test_id': test_id,
            'status': 'error',
            'error': str(e)
        })
        raise

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

                while not stop_event.is_alive():
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
            # Check if database exists and print diagnostic info
            print(f"\nChecking database at: {self.db_path}")
            
            # Detailed file check
            try:
                file_exists = os.path.exists(self.db_path)
                file_size = os.path.getsize(self.db_path) if file_exists else 0
                file_permissions = oct(os.stat(self.db_path).st_mode)[-3:] if file_exists else "N/A"
                print(f"File exists: {file_exists}")
                print(f"File size: {file_size} bytes")
                print(f"File permissions: {file_permissions}")
            except Exception as e:
                print(f"Error checking file: {e}")

            if not os.path.exists(self.db_path):
                print("\nNo rankings available - database file not found.")
                print(f"Expected database at: {self.db_path}")
                print("Please run some tests first using --test flag to create the database.")
                logging.warning(f"Database file not found at: {self.db_path}")
                return
            
            try:
                with sqlite3.connect(self.db_path) as conn:
                    # First check if the required tables exist
                    c = conn.cursor()
                    
                    # Get all tables in the database
                    c.execute("SELECT name FROM sqlite_master WHERE type='table'")
                    all_tables = c.fetchall()
                    print("\nAll tables in database:")
                    for table in all_tables:
                        print(f"- {table[0]}")
                        # Count rows in each table
                        try:
                            c.execute(f"SELECT COUNT(*) FROM {table[0]}")
                            count = c.fetchone()[0]
                            print(f"  Rows: {count}")
                        except sqlite3.Error as e:
                            print(f"  Error counting rows: {e}")
                    
                    # Check for our specific tables
                    c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name IN ('solver_versions', 'test_results')")
                    tables = c.fetchall()
                    tables = [t[0] for t in tables]
                    
                    if 'solver_versions' not in tables or 'test_results' not in tables:
                        print("\nDatabase exists but required tables are missing.")
                        print(f"Found tables: {', '.join(tables)}")
                        print("Please run some tests first using --test flag to initialize the database.")
                        return

                    # Check solver_versions content
                    print("\nSolver versions in database:")
                    c.execute("SELECT version_id, version_name, description FROM solver_versions")
                    versions = c.fetchall()
                    for v in versions:
                        print(f"- ID: {v[0]}, Name: {v[1]}, Desc: {v[2]}")

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

                    if not results:
                        print("\nNo rankings available - no test results found in database.")
                        print("Please run some tests first using --test flag.")
                        logging.warning("No test results found in database")
                        return

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
                    
            except sqlite3.OperationalError as e:
                print(f"\nError accessing database: {e}")
                print("The database file might be corrupted or not properly initialized.")
                print("Try running some tests first using --test flag.")
                logging.error(f"SQLite operational error: {e}")
                
        except Exception as e:
            error_msg = f"Error displaying rankings: {e}"
            print(error_msg)
            logging.error(error_msg)
            if args.debug:
                import traceback
                traceback.print_exc()

    def debug_processes(self, solver_path: str, num_processes: int = 1) -> None:
        """
        Debug mode to track process creation, initialization and termination.
        """
        logging.info("Starting process debugging mode")
        self.solver_path = solver_path
        
        if not os.path.exists(solver_path):
            logging.error(f"Solver not found at path: {solver_path}")
            return

        # Create a queue for inter-process communication
        result_queue = multiprocessing.Queue()
        processes = []
        
        try:
            # Initialize display manager for debug output
            if self.use_ui:
                # For debug mode, we'll use 1 test per process
                total_tests = num_processes
                self.display_manager = UnifiedDisplayManager(total_tests, num_processes, use_ui=self.use_ui)
                logging.info("Display manager initialized")
            
            # Create test instance for debugging
            test_instance = self.test_dir / "easy/test_1/instance.txt"
            if not test_instance.exists():
                logging.error(f"Test instance not found at: {test_instance}")
                return
            
            logging.info(f"Using test instance: {test_instance}")
            
            # Start processes one by one with detailed logging
            for i in range(num_processes):
                process_id = i + 1
                logging.info(f"Initializing process {process_id}")
                
                try:
                    p = multiprocessing.Process(
                        target=debug_subprocess,
                        args=(
                            solver_path,
                            test_instance,
                            1,  # version_id
                            1,  # test_id
                            process_id,
                            self.db_path,
                            result_queue
                        )
                    )
                    processes.append(p)
                    logging.info(f"Created process object {process_id}")
                    
                    p.start()
                    logging.info(f"Started process {process_id} with PID {p.pid}")
                    
                    if self.display_manager:
                        self.display_manager.update_process(
                            process_id - 1,  # UI uses 0-based indexing
                            0.0,  # initial progress
                            "Starting",
                            0.0,  # initial fitness
                            0.0,  # initial coverage
                            0.0,  # initial edge score
                            0.0   # initial theoretical max
                        )
                        self.display_manager.update_test_info(process_id - 1, str(test_instance))
                    
                    # Wait a bit between process starts to avoid race conditions
                    time.sleep(0.5)
                    
                except Exception as e:
                    logging.error(f"Failed to create process {process_id}: {str(e)}")
                    raise
            
            # Monitor process states
            while processes:
                for p in processes[:]:
                    if not p.is_alive():
                        logging.info(f"Process {p.pid} has terminated")
                        p.join()
                        processes.remove(p)
                        
                # Check for results
                while not result_queue.empty():
                    result = result_queue.get()
                    process_id = result['process_id']
                    logging.info(f"Received result from process {process_id}: {result['status']}")
                    
                    if self.display_manager:
                        status = "Completed" if result['status'] == 'success' else "Error"
                        progress = 100.0 if result['status'] == 'success' else 0.0
                        fitness = float(result.get('result', 0.0))
                        
                        self.display_manager.update_process(
                            process_id - 1,  # UI uses 0-based indexing
                            progress,
                            status,
                            fitness
                        )
                        self.display_manager.increment_completed()
                
                time.sleep(0.1)
                
        except KeyboardInterrupt:
            logging.info("Received keyboard interrupt, cleaning up processes...")
            for p in processes:
                if p.is_alive():
                    logging.info(f"Terminating process {p.pid}")
                    p.terminate()
                    p.join()
        
        finally:
            if self.display_manager:
                self.display_manager.cleanup()
                logging.info("Display manager cleaned up")
            logging.info("Process debugging completed")

def debug_subprocess(
        solver_path: str,
        instance_path: Path,
        version_id: int,
        test_id: int,
        process_id: int,
        db_path: str,
        queue: multiprocessing.Queue
) -> None:
    """
    Debug version of solver subprocess with enhanced error handling and resource management
    """
    logger = setup_logger(True)
    
    try:
        logger.info(f"Process {process_id}: Debug subprocess started")
        
        # Convert to absolute paths and validate files
        solver_abs_path = os.path.abspath(solver_path)
        instance_abs_path = os.path.abspath(str(instance_path))
        
        # Use config from build directory
        build_dir = os.path.dirname(solver_abs_path)
        config_file = os.path.join(build_dir, "config.cfg")
        
        # Validate solver executable
        if not os.path.exists(solver_abs_path):
            raise FileNotFoundError(f"Solver not found: {solver_abs_path}")
        
        # Check if solver is executable
        if not os.access(solver_abs_path, os.X_OK):
            logger.error(f"Process {process_id}: Solver is not executable")
            os.chmod(solver_abs_path, 0o755)  # Try to make it executable
            if not os.access(solver_abs_path, os.X_OK):
                raise PermissionError(f"Cannot make solver executable: {solver_abs_path}")
        
        # Validate instance file
        if not os.path.exists(instance_abs_path):
            raise FileNotFoundError(f"Instance file not found: {instance_abs_path}")
            
        # Validate config file
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file not found: {config_file}")
            
        # Create process directory
        process_dir = Path(f"process_{process_id}")
        process_dir.mkdir(exist_ok=True)
        
        # Prepare output file
        output_file = process_dir / f"debug_result_{process_id}.txt"
        
        # Log all paths and command that will be executed
        logger.info(f"Process {process_id}: Using paths:")
        logger.info(f"- Solver: {solver_abs_path}")
        logger.info(f"- Instance: {instance_abs_path}")
        logger.info(f"- Config: {config_file}")
        logger.info(f"- Output: {output_file}")
        
        # Prepare command with all required arguments
        cmd = [
            solver_abs_path,
            "test_instance",
            "-i", str(instance_abs_path),
            "-o", str(output_file),
            "-cfg", str(config_file),
            "-pid", str(process_id),
            "-diff", "Easy"
        ]
        
        logger.info(f"Process {process_id}: Executing command: {' '.join(cmd)}")
        
        # Start solver process with enhanced error handling
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            bufsize=1
        )
        
        # Monitor process output in real-time
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                line = output.strip()
                logger.debug(f"Process {process_id} stdout: {line}")
                
                # Check for specific error indicators in output
                if "Failed to load instance" in line:
                    logger.error(f"Process {process_id}: Instance loading failed")
                elif "Invalid format" in line:
                    logger.error(f"Process {process_id}: Invalid instance format")
                elif "Memory allocation failed" in line:
                    logger.error(f"Process {process_id}: Memory allocation error")
                    
        # Get return code and remaining output
        return_code = process.wait()
        _, stderr = process.communicate()
        
        if stderr:
            logger.error(f"Process {process_id} stderr: {stderr}")
            
        if return_code < 0:
            logger.error(f"Process {process_id}: Process terminated with signal {-return_code}")
            raise RuntimeError(f"Process terminated with signal {-return_code}")
        elif return_code > 0:
            logger.error(f"Process {process_id}: Process failed with return code {return_code}")
            raise RuntimeError(f"Process failed with return code {return_code}")
            
        logger.info(f"Process {process_id}: Process completed successfully")
        
    except Exception as e:
        logger.error(f"Process {process_id}: Fatal error: {str(e)}")
        queue.put({
            'process_id': process_id,
            'version_id': version_id,
            'test_id': test_id,
            'status': 'error',
            'error': str(e)
        })
        raise
    finally:
        try:
            if 'process_dir' in locals() and process_dir.exists():
                shutil.rmtree(process_dir)
                logger.info(f"Process {process_id}: Cleaned up process directory")
        except Exception as e:
            logger.error(f"Process {process_id}: Cleanup error: {str(e)}")

def monitor_solver_process(
    solver: SolverProcess,
    queue: multiprocessing.Queue,
    process_id: int,
    version_id: int,
    test_id: int
) -> None:
    """Monitor solver process with timeout handling"""
    logger = logging.getLogger()
    stdout_buffer = []
    stderr_buffer = []
    fitness = 0.0
    
    while True:
        # Check for timeout
        if solver.check_timeout():
            raise ProcessTimeoutError(f"Process {process_id} exceeded timeout of {PROCESS_TIMEOUT} seconds")
            
        # Read output with timeout
        try:
            stdout = solver.process.stdout.readline()
            stderr = solver.process.stderr.readline()
            
            if stdout == '' and stderr == '' and solver.process.poll() is not None:
                break
                
            if stdout:
                handle_stdout_line(stdout.strip(), stdout_buffer, process_id, fitness, queue, version_id, test_id)
                
            if stderr:
                handle_stderr_line(stderr.strip(), stderr_buffer, process_id)
                
        except Exception as e:
            logger.error(f"Process {process_id}: Error reading output: {e}")
            break
            
    # Process completion
    return_code = solver.process.wait(timeout=1.0)
    if return_code != 0:
        raise RuntimeError(f"Process {process_id} failed with return code {return_code}")

def handle_stdout_line(
    line: str,
    buffer: List[str],
    process_id: int,
    fitness: float,
    queue: multiprocessing.Queue,
    version_id: int,
    test_id: int
) -> None:
    """Handle a line of stdout from the solver process"""
    logger = logging.getLogger()
    buffer.append(line)
    logger.debug(f"Process {process_id} stdout: {line}")
    
    # Parse fitness if present
    if "Best Fitness:" in line:
        try:
            fitness = float(line.split("Best Fitness:")[1].strip())
            logger.info(f"Process {process_id}: Current best fitness: {fitness}")
            
            # Send progress update
            queue.put({
                'process_id': process_id,
                'version_id': version_id,
                'test_id': test_id,
                'status': 'running',
                'result': fitness
            })
        except (ValueError, IndexError):
            logger.warning(f"Process {process_id}: Could not parse fitness from line: {line}")

def handle_stderr_line(line: str, buffer: List[str], process_id: int) -> None:
    """Handle a line of stderr from the solver process"""
    logger = logging.getLogger()
    buffer.append(line)
    logger.error(f"Process {process_id} stderr: {line}")

def main() -> None:
    parser = argparse.ArgumentParser(description="SBH Testing Framework (multiprocessing + curses)")

    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define paths relative to the parent directory of scripts
    parent_dir = os.path.dirname(script_dir)
    default_solver = os.path.join(parent_dir, 'optymalizacja_kombinatoryczna')
    default_db = os.path.join(parent_dir, 'results.db')
    default_test_dir = os.path.join(parent_dir, 'test_cases')

    parser.add_argument('--generate', action='store_true', help='Generate test cases (if not exist)')
    parser.add_argument('--test', action='store_true', help='Run tests on solver in parallel')
    parser.add_argument('--solver', default=str(default_solver), help='Path to solver executable')
    parser.add_argument('--version', help='Solver version name')
    parser.add_argument('--desc', help='Solver version description')
    parser.add_argument('--processes', type=int, default=1, help='Number of parallel processes')
    parser.add_argument('--rankings', action='store_true', help='Display solver rankings')
    parser.add_argument('--log', action='store_true', help='Use logging instead of UI')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--debug-processes', action='store_true', help='Run in process debugging mode')
    parser.add_argument('--db', default=str(default_db), help='Path to results database')

    args = parser.parse_args()

    try:
        use_ui = not args.log and sys.stdout.isatty()
        logger = setup_logger(debug_mode=args.debug)
        
        # Print diagnostic information
        print(f"Current working directory: {os.getcwd()}")
        print(f"Script directory: {script_dir}")
        print(f"Database path: {args.db}")
        print(f"Database exists: {os.path.exists(args.db)}")
        
        if args.debug:
            logging.info(f"Current working directory: {os.getcwd()}")
            logging.info(f"Script directory: {script_dir}")
            logging.info(f"Database path: {args.db}")
            logging.info(f"Database exists: {os.path.exists(args.db)}")

        framework = SBHTestingFramework(
            db_path=args.db,
            test_dir=str(default_test_dir),
            use_ui=use_ui
        )
        global framework_instance
        framework_instance = framework

        # Only check for solver executable if we need it
        if args.test or args.debug_processes or args.generate:
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

        if args.debug_processes:
            logging.info(f"Starting debug mode with {args.processes} processes")
            framework.debug_processes(args.solver, args.processes)
        elif args.test:
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
