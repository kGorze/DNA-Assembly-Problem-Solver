# display_manager.py

import curses
import threading
import logging
import time
import sys
import signal
from typing import Optional, Dict, Any


class UnifiedDisplayManager:
    def __init__(self, total_tests: int, num_processes: int, use_ui: bool = True) -> None:
        """
        Initialize the UnifiedDisplayManager.

        Args:
            total_tests (int): Total number of test cases.
            num_processes (int): Number of parallel processes.
            use_ui (bool): Whether to use the UI (curses) or not.
        """
        self.total_tests = total_tests
        self.num_processes = num_processes
        self.use_ui = use_ui
        self.output_lock = threading.Lock()
        self.output_messages = []
        self.process_status: Dict[int, Dict[str, Any]] = {
            pid: {"progress": 0.0, "status": "Initializing", "fitness": 0.0} for pid in range(num_processes)
        }
        self.completed_tests = 0
        self.screen: Optional[curses.window] = None
        self.window_lock = threading.Lock()

        if self.use_ui:
            # Start the curses UI in a separate thread
            self.ui_thread = threading.Thread(target=self._start_curses_ui, daemon=True)
            self.ui_thread.start()
            logging.info("UnifiedDisplayManager: Curses UI thread started.")
        else:
            logging.info("UnifiedDisplayManager: UI disabled. Logging only.")

    def _start_curses_ui(self) -> None:
        """Initialize and manage the curses UI."""
        try:
            self.screen = curses.initscr()
            curses.noecho()
            curses.cbreak()
            self.screen.nodelay(True)  # Make getch() non-blocking
            curses.start_color()
            curses.init_pair(1, curses.COLOR_RED, curses.COLOR_BLACK)    # For errors
            curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)  # For info

            self._draw_initial_screen()

            while True:
                with self.window_lock:
                    self._update_process_windows()
                    self._update_output_window()
                    self._update_footer()
                time.sleep(0.1)  # Refresh rate

        except Exception as e:
            logging.exception(f"UnifiedDisplayManager: Exception in curses UI: {e}")
        finally:
            if self.use_ui and self.screen:
                curses.nocbreak()
                self.screen.keypad(False)
                curses.echo()
                curses.endwin()
                logging.info("UnifiedDisplayManager: Curses UI cleaned up.")

    def _draw_initial_screen(self) -> None:
        """Draw the initial layout of the UI."""
        self.screen.clear()
        height, width = self.screen.getmaxyx()

        # Divide the screen into sections
        process_win_height = 3
        output_win_height = height - (process_win_height * self.num_processes) - 4
        process_windows = []
        for pid in range(self.num_processes):
            win = curses.newwin(process_win_height, width, pid * process_win_height, 0)
            win.box()
            win.addstr(0, 2, f" Process {pid} ")
            process_windows.append(win)
            win.refresh()
            logging.debug(f"UnifiedDisplayManager: Process window {pid} created.")
        self.process_windows = process_windows

        # Output window
        self.output_window = curses.newwin(output_win_height, width, self.num_processes * process_win_height, 0)
        self.output_window.box()
        self.output_window.addstr(0, 2, " Output ")
        self.output_window.scrollok(True)
        self.output_window.refresh()
        logging.debug("UnifiedDisplayManager: Output window created.")

        # Footer
        self.footer_window = curses.newwin(2, width, self.num_processes * process_win_height + output_win_height, 0)
        self.footer_window.box()
        self.footer_window.addstr(0, 2, " Footer ")
        self.footer_window.addstr(1, 2, f"Completed Tests: {self.completed_tests}/{self.total_tests}")
        self.footer_window.refresh()
        logging.debug("UnifiedDisplayManager: Footer window created.")

    def _update_process_windows(self) -> None:
        """Update the status of each process in the UI."""
        for pid, status in self.process_status.items():
            win = self.process_windows[pid]
            win.erase()
            win.box()
            win.addstr(0, 2, f" Process {pid} ")

            progress_str = f"Progress: {status['progress']:.2f}%"
            status_str = f"Status: {status['status']}"
            fitness_str = f"Fitness: {status['fitness']:.2f}"

            win.addstr(1, 2, progress_str)
            win.addstr(1, 25, status_str)
            win.addstr(1, 50, fitness_str)
            win.refresh()
            logging.debug(f"UnifiedDisplayManager: Updated process {pid} window.")

    def _update_output_window(self) -> None:
        """Update the output window with new messages."""
        with self.output_lock:
            for message in self.output_messages:
                if self.use_ui:
                    self.output_window.addstr(message + "\n")
                logging.debug(f"UnifiedDisplayManager: Added output message: {message}")
            self.output_messages.clear()
            if self.use_ui:
                self.output_window.refresh()

    def _update_footer(self) -> None:
        """Update the footer with the number of completed tests."""
        if self.use_ui:
            self.footer_window.erase()
            self.footer_window.box()
            self.footer_window.addstr(0, 2, " Footer ")
            self.footer_window.addstr(1, 2, f"Completed Tests: {self.completed_tests}/{self.total_tests}")
            self.footer_window.refresh()
            logging.debug("UnifiedDisplayManager: Updated footer window.")

    def add_output(self, message: str, error: bool = False) -> None:
        """
        Add a message to the output window.

        Args:
            message (str): The message to display.
            error (bool): Whether the message is an error.
        """
        with self.output_lock:
            self.output_messages.append(message)
            if self.use_ui:
                # Optionally, color the message based on error flag
                try:
                    if error:
                        self.output_window.attron(curses.color_pair(1))
                        self.output_window.addstr(message + "\n")
                        self.output_window.attroff(curses.color_pair(1))
                    else:
                        self.output_window.attron(curses.color_pair(2))
                        self.output_window.addstr(message + "\n")
                        self.output_window.attroff(curses.color_pair(2))
                except curses.error as e:
                    logging.error(f"UnifiedDisplayManager: Curses error while adding output: {e}")
            logging.info(f"UnifiedDisplayManager: add_output called with message: {message}, error: {error}")

    def update_process(self, pid: int, progress: float, status: str, fitness: float) -> None:
        """
        Update the status of a specific process.

        Args:
            pid (int): Process ID.
            progress (float): Progress percentage.
            status (str): Current status message.
            fitness (float): Fitness score.
        """
        if pid not in self.process_status:
            logging.warning(f"UnifiedDisplayManager: Received update for unknown pid {pid}.")
            return
        self.process_status[pid]['progress'] = progress
        self.process_status[pid]['status'] = status
        self.process_status[pid]['fitness'] = fitness
        logging.info(f"UnifiedDisplayManager: Process {pid} updated - Progress: {progress}%, Status: {status}, Fitness: {fitness}")

    def increment_completed(self) -> None:
        """
        Increment the count of completed tests.
        """
        self.completed_tests += 1
        logging.info(f"UnifiedDisplayManager: Incremented completed_tests to {self.completed_tests}/{self.total_tests}")

    def cleanup(self) -> None:
        if self.use_ui:
            try:
                # Najpierw wyłącz flagę UI
                self.use_ui = False
                # Poczekaj na zakończenie wątku UI
                time.sleep(0.2)

                # Przywróć standardowe ustawienia terminala
                if self.screen:
                    self.screen.keypad(False)
                    curses.nocbreak()
                    curses.echo()
                    curses.endwin()

                logging.info("UnifiedDisplayManager: Cleanup completed successfully")
            except Exception as e:
                logging.exception(f"UnifiedDisplayManager: Exception during cleanup: {e}")
    
    def wait_for_ui(self) -> None:
        """
        Wait for the UI thread to finish (useful for non-daemon threads).
        """
        if self.use_ui and hasattr(self, 'ui_thread'):
            self.ui_thread.join()
            logging.info("UnifiedDisplayManager: UI thread has been joined.")
