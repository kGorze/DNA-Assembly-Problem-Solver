# sbh_testing_framework.py

#!/usr/bin/env python3

"""
sbh_testing_framework.py - Enhanced Testing Framework for SBH Solvers with CLI Support
using multiprocessing.Queue to communicate between subprocesses and the main process for real-time curses updates.
"""

import os
import sys
import sqlite3
import time
import argparse
import subprocess
import multiprocessing
import threading
import logging
import signal
import re
import curses
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any
from display_manager import UnifiedDisplayManager


# Global reference to the framework instance, used by the signal handler
framework_instance = None

# Constants for default paths and configurations
DEFAULT_LOG_FILENAME = "sbh_framework.log"
DEFAULT_BUILD_DIR = 'build'
DEFAULT_SOLVER = os.path.join(DEFAULT_BUILD_DIR, 'optymalizacja_kombinatoryczna')
DEFAULT_DB = os.path.join(DEFAULT_BUILD_DIR, 'results.db')
DEFAULT_TEST_DIR = os.path.join(DEFAULT_BUILD_DIR, 'test_cases')
MAX_LOG_SIZE = 5 * 1024 * 1024  # 5 MB
BACKUP_LOG_COUNT = 5

# Configure logger with file rotation
from logging.handlers import RotatingFileHandler

handler = RotatingFileHandler(DEFAULT_LOG_FILENAME, maxBytes=MAX_LOG_SIZE, backupCount=BACKUP_LOG_COUNT)
logging.basicConfig(
    level=logging.DEBUG,  # Set to DEBUG for detailed logs
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        handler
    ]
)

#handlers=[
#handler,
#logging.StreamHandler()
#]

# Test log message
logging.info("Logging setup complete. Logs will be saved to: {}".format(os.path.abspath(DEFAULT_LOG_FILENAME)))

def cleanup_and_exit(signum, frame):
    if framework_instance and framework_instance.display_manager:
        framework_instance.display_manager.cleanup()
    # Przywróć standardowe ustawienia terminala
    try:
        curses.endwin()
    except:
        pass
    sys.exit(0)

signal.signal(signal.SIGINT, cleanup_and_exit)

# Make sure framework_instance is declared as global at the start of the file
framework_instance = None

def signal_handler(signum: int, frame: Any) -> None:
    """Handle interruption signals gracefully."""
    global framework_instance
    logging.info(f"Received signal {signum}. Initiating cleanup.")
    if framework_instance and framework_instance.display_manager:
        framework_instance.display_manager.cleanup()
    print("\nReceived interrupt signal. Cleaning up...")
    sys.exit(1)

signal.signal(signal.SIGINT, signal_handler)

def solver_subprocess(
        solver_path: str,
        instance_path: Path,
        version_id: int,
        test_id: int,
        process_id: int,
        db_path: str,
        queue: multiprocessing.Queue
) -> None:
    """
    Function to run in each subprocess:
    1. Start the solver.
    2. Read stdout line by line.
       - If line starts with "PROGRESS_UPDATE:", parse and send to queue.
       - Else, send as output line.
    3. After solver finishes, read results.txt and send 'done' message.
    
    Args:
        solver_path (str): Path to the solver executable.
        instance_path (Path): Path to the instance file.
        version_id (int): ID of the solver version.
        test_id (int): ID of the test case.
        process_id (int): ID of the process.
        db_path (str): Path to the database file.
        queue (multiprocessing.Queue): Queue for inter-process communication.
    """
    start_time = time.time()
    logging.info(f"Process {process_id}: Starting solver for test_id {test_id}.")

    # Regular expression to parse progress lines
    progress_pattern = re.compile(r"PROGRESS_UPDATE:(\d+):([0-9\.]+):(.+):([0-9\.]+)$")

    try:
        cmd = [
            solver_path,
            "test_instance",
            "-i", str(instance_path),
            "-o", str(instance_path.with_name("results.txt")),
            "-pid", str(process_id)
        ]

        logging.info(f"Process {process_id}: Executing command: {' '.join(cmd)}")
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1
        )

        # Notify that the solver has started
        queue.put(("started", process_id))
        logging.info(f"Process {process_id}: Solver started.")

        # Read stdout line by line
        for line in proc.stdout:
            line = line.strip()
            if line.startswith("PROGRESS_UPDATE:"):
                # Use regex to parse the progress line
                match = progress_pattern.match(line)
                if match:
                    try:
                        pid = int(match.group(1))
                        progress = float(match.group(2))
                        status = match.group(3)
                        fitness = float(match.group(4))
                        # Send progress update to main process
                        queue.put(("progress", pid, progress, status, fitness))
                        logging.debug(f"Process {process_id}: Progress update - PID: {pid}, Progress: {progress}, Status: {status}, Fitness: {fitness}")
                    except (ValueError, IndexError) as e:
                        error_msg = f"Process {process_id}: Error converting progress update: {e}"
                        queue.put(("output", process_id, error_msg, True))
                        logging.error(error_msg)
                else:
                    error_msg = f"Process {process_id}: Unrecognized progress update format: {line}"
                    queue.put(("output", process_id, error_msg, True))
                    logging.error(error_msg)
            else:
                # Regular output line
                queue.put(("output", process_id, line, False))
                logging.debug(f"Process {process_id}: Output - {line}")

        proc.wait()
        exit_code = proc.returncode
        logging.info(f"Process {process_id}: Solver exited with code {exit_code}.")

        if exit_code != 0:
            error_msg = f"Process {process_id}: Solver exited with code {exit_code}"
            queue.put(("output", process_id, error_msg, True))
            logging.error(error_msg)

        # Read results.txt if it exists
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
                        logging.info(f"Process {process_id}: Results - Original DNA: {original_dna}, Reconstructed DNA: {reconstructed_dna}, Distance: {dist}")
            else:
                logging.warning(f"Process {process_id}: Results file not found: {results_path}")
        except Exception as e:
            error_msg = f"Process {process_id}: Error reading results file: {e}"
            queue.put(("output", process_id, error_msg, True))
            logging.error(error_msg)

        exec_time_ms = int((time.time() - start_time) * 1000)
        logging.info(f"Process {process_id}: Execution time: {exec_time_ms} ms.")

        # Send 'done' message to main process
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
        """
        Initialize the SBH Testing Framework with improved error handling.

        Args:
            db_path (str): Path to the SQLite database.
            test_dir (str): Directory where test cases are stored.
            use_ui (bool): Whether to use the UI (curses) or not.
        """
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
        """Initialize SQLite database with required schema."""
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
        """
        Generate test cases using the C++ program's generate_instance mode.

        It generates test cases for 'easy', 'medium', and 'hard' difficulties with specified parameters.
        """
        if self.test_dir.exists() and any(self.test_dir.iterdir()):
            print("Test cases already exist, skipping generation.")
            logging.info("Test cases already exist, skipping generation.")
            return

        os.makedirs(self.test_dir, exist_ok=True)
        logging.info(f"Created test directory at {self.test_dir}.")

        # Define test configurations for different difficulty levels
        test_configs: Dict[str, Dict[str, Any]] = {
            "easy": {
                "counts": 20,
                "params": [
                    {"n": 300, "k": 7, "dk": 0, "ln": 10, "lp": 10},
                    {"n": 325, "k": 7, "dk": 0, "ln": 11, "lp": 11},
                    {"n": 350, "k": 8, "dk": 0, "ln": 12, "lp": 12},
                    {"n": 375, "k": 8, "dk": 0, "ln": 13, "lp": 13}
                ]
            },
            "medium": {
                "counts": 15,
                "params": [
                    {"n": 400, "k": 8, "dk": 1, "ln": 15, "lp": 15},
                    {"n": 425, "k": 8, "dk": 1, "ln": 16, "lp": 16},
                    {"n": 450, "k": 9, "dk": 1, "ln": 18, "lp": 18},
                    {"n": 475, "k": 9, "dk": 1, "ln": 19, "lp": 19}
                ]
            },
            "hard": {
                "counts": 10,
                "params": [
                    {"n": 500, "k": 9, "dk": 2, "ln": 20, "lp": 20},
                    {"n": 550, "k": 9, "dk": 2, "ln": 22, "lp": 22},
                    {"n": 600, "k": 10, "dk": 2, "ln": 25, "lp": 25},
                    {"n": 650, "k": 10, "dk": 2, "ln": 27, "lp": 27}
                ]
            }
        }

        # Validate solver path existence
        if not self.solver_path or not os.path.isfile(self.solver_path):
            error_msg = f"Solver executable not found at path: {self.solver_path}"
            print(error_msg)
            logging.error(error_msg)
            sys.exit(1)
        logging.info(f"Validated solver path: {self.solver_path}")

        # Generate test cases using the C++ program
        with sqlite3.connect(self.db_path) as conn:
            for difficulty, config in test_configs.items():
                diff_dir = self.test_dir / difficulty
                diff_dir.mkdir(parents=True, exist_ok=True)
                logging.info(f"Created directory for difficulty '{difficulty}' at {diff_dir}.")

                for i in range(config["counts"]):
                    for params in config["params"]:
                        test_subdir = diff_dir / f"test_{i + 1}"
                        test_subdir.mkdir(parents=True, exist_ok=True)

                        instance_path = test_subdir / "instance.txt"

                        # Generate instance using the C++ program
                        cmd = (
                            f"{self.solver_path} generate_instance "
                            f"-n {params['n']} -k {params['k']} -dk {params['dk']} "
                            f"-ln {params['ln']} -lp {params['lp']} "
                            f"-o {instance_path}"
                        )

                        try:
                            logging.debug(f"Executing command: {cmd}")
                            subprocess.run(cmd, shell=True, check=True)
                            logging.info(f"Generated instance at {instance_path}.")
                        except subprocess.CalledProcessError as e:
                            error_msg = f"Error generating instance {instance_path}: {e}"
                            print(error_msg)
                            logging.error(error_msg)
                            continue

                        # Save test case information to the database
                        try:
                            conn.execute("""
                                INSERT INTO test_cases (difficulty, instance_path, created_at)
                                VALUES (?, ?, ?)
                            """, (difficulty, str(instance_path), datetime.now()))
                            logging.debug(f"Inserted test case into database: {instance_path}")
                        except sqlite3.Error as e:
                            error_msg = f"Database error while inserting test case {instance_path}: {e}"
                            print(error_msg)
                            logging.error(error_msg)

            conn.commit()
            logging.info("All test cases generated and saved to the database.")

    def test_solver_version(
            self,
            solver_path: str,
            version_name: str,
            description: str,
            num_processes: Optional[int] = None
    ) -> None:
            """
            Run solver in multiple subprocesses in parallel,
            using multiprocessing.Queue for real-time communication with main process.
        
            Args:
                solver_path (str): Path to the solver executable.
                version_name (str): Name of the solver version.
                description (str): Description of the solver version.
                num_processes (Optional[int]): Number of parallel processes.
            """
            logging.info(f"Starting tests for solver version: {version_name}")
            try:
                # Register the solver version in the database
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
        
                # Get list of test_cases from the database
                with sqlite3.connect(self.db_path) as conn:
                    c = conn.cursor()
                    c.execute("SELECT test_id, instance_path FROM test_cases")
                    test_cases = c.fetchall()
        
                if not test_cases:
                    message = "No test cases found. Please generate first."
                    if self.use_ui and self.display_manager:
                        self.display_manager.add_output(message, error=True)
                    else:
                        logging.error(message)
                    return
        
                total_tests = len(test_cases)
                logging.info(f"Total test cases to run: {total_tests}")
        
                # Determine the number of concurrent processes
                cpu_count = os.cpu_count() or 1
                if num_processes is None:
                    num_processes = cpu_count
                else:
                    num_processes = min(num_processes, cpu_count, total_tests)
                logging.info(f"Using {num_processes} parallel processes out of {cpu_count} available CPU cores.")
        
                # Initialize the display manager
                self.display_manager = UnifiedDisplayManager(
                    total_tests=total_tests,
                    num_processes=num_processes,
                    use_ui=self.use_ui
                )
                logging.info("Initialized UnifiedDisplayManager.")
        
                # Extend the add_output method to include logging
                if self.display_manager:
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
        
                # Create a multiprocessing Queue for communication
                manager = multiprocessing.Manager()
                queue = manager.Queue()
                logging.info("Created multiprocessing.Queue for inter-process communication.")
        
                # Prepare the list of tests and initialize process tracking
                tests_to_run = list(test_cases)
                running_processes: list[multiprocessing.Process] = []
                test_index = 0
                logging.info("Prepared list of tests to run.")
        
                # Start the initial batch of processes
                for pid in range(num_processes):
                    if test_index >= total_tests:
                        break
                    test_id, inst_path = tests_to_run[test_index]
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
                    logging.info(f"Started process {pid} for test_id {test_id}.")
                    test_index += 1
        
                # Start a thread to read from the queue
                stop_event = threading.Event()
        
                def queue_listener() -> None:
                    nonlocal running_processes, test_index, total_tests, tests_to_run
                    finished_tests = 0
        
                    while not stop_event.is_set() and finished_tests < total_tests:
                        try:
                            msg = queue.get(timeout=1)
                        except Exception as e:
                            logging.error(f"Queue listener: Exception while getting message: {e}")
                            continue
        
                        if not msg:
                            continue
        
                        mtype = msg[0]
        
                        if mtype == "progress":
                            _, pid, progress, status, fitness = msg
                            self.display_manager.update_process(pid, progress, status, fitness)
                            logging.debug(f"Main process: Received progress from PID {pid}: {progress}% - {status} - Fitness: {fitness}")
        
                        elif mtype == "output":
                            _, pid, line, is_error = msg
                            self.display_manager.add_output(line, error=is_error)
                            if is_error:
                                logging.error(f"Process {pid} Output: {line}")
                            else:
                                logging.info(f"Process {pid} Output: {line}")
        
                        elif mtype == "done":
                            _, pid, test_id, dist, orig, recon, exec_ms = msg
                            logging.info(f"Main process: Received 'done' from PID {pid} for test_id {test_id} with distance {dist} in {exec_ms} ms.")
        
                            try:
                                # Save to database
                                with sqlite3.connect(self.db_path) as conn:
                                    c = conn.cursor()
                                    c.execute("""
                                            INSERT INTO test_results
                                            (version_id, test_id, levenshtein_distance, 
                                             execution_time_ms, original_dna, reconstructed_dna, tested_at)
                                            VALUES (?, ?, ?, ?, ?, ?, ?)
                                        """, (
                                        version_id,
                                        test_id,
                                        dist,
                                        exec_ms,
                                        orig,
                                        recon,
                                        datetime.now()
                                    ))
                                    conn.commit()
                                logging.info(f"Main process: Saved results for test_id {test_id} to database.")
                            except sqlite3.Error as e:
                                error_msg = f"Main process: Database error while saving results for test_id {test_id}: {e}"
                                self.display_manager.add_output(error_msg, error=True)
                                logging.error(error_msg)
        
                            self.display_manager.add_output(
                                f"Test {test_id} (process {pid}) finished. dist={dist}",
                                error=False
                            )
                            self.display_manager.increment_completed()
                            finished_tests += 1
        
                            logging.info(f"Main process: Finished_tests count: {finished_tests}/{total_tests}")
        
                            # Start a new process if there are remaining tests
                            if test_index < total_tests:
                                new_test_id, new_inst_path = tests_to_run[test_index]
                                p = multiprocessing.Process(
                                    target=solver_subprocess,
                                    args=(
                                        solver_path,
                                        Path(new_inst_path),
                                        version_id,
                                        new_test_id,
                                        pid,  # Reuse the same process_id
                                        self.db_path,
                                        queue
                                    )
                                )
                                p.start()
                                running_processes.append(p)
                                logging.info(f"Started new process {pid} for test_id {new_test_id}.")
                                test_index += 1
        
                            # Log the current number of running processes
                            active_procs = [p for p in running_processes if p.is_alive()]
                            logging.info(f"Main process: Currently running processes: {len(active_procs)}")
        
                    # After finishing the loop
                    logging.info("Queue listener has finished processing all tests.")
        
                # Start the queue listener thread
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
                    raise  # Re-raise to be caught by outer try block
            except Exception as e:
                logging.error(f"Error in test_solver_version: {e}")
                raise  # Re-raise after logging
            finally:
                # Cleanup code all in one finally block
                stop_event.set()
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
        """
        Display solver rankings with improved formatting.
        """
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
    """Main entry point with improved error handling."""
    parser = argparse.ArgumentParser(description="SBH Testing Framework (multiprocessing + curses)")

    script_dir = Path(__file__).parent
    build_dir = script_dir / DEFAULT_BUILD_DIR
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

    args = parser.parse_args()

    try:
        # Determine if UI should be used
        use_ui = not args.log and sys.stdout.isatty()
        logging.info(f"Determined use_ui: {use_ui}")

        # Create framework instance
        framework = SBHTestingFramework(
            db_path=str(default_db),
            test_dir=str(default_test_dir),
            use_ui=use_ui
        )
        # Set global reference to framework for signal handler
        global framework_instance
        framework_instance = framework

        framework.solver_path = args.solver
        logging.info(f"Using solver at path: {args.solver}")

        # Validate solver path existence
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
