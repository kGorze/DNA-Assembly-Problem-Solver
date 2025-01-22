import curses
import threading
import logging
import time
import sys
from typing import Optional, Dict, Any
from datetime import datetime, timedelta

class UnifiedDisplayManager:
    def __init__(self, total_tests: int, num_processes: int, use_ui: bool = True) -> None:
        self.total_tests = total_tests
        self.num_processes = num_processes
        self.use_ui = use_ui
        self.output_lock = threading.Lock()
        self.output_messages = []
        self.process_status: Dict[int, Dict[str, Any]] = {
            pid: {
                "progress": 0.0,
                "status": "Initializing",
                "fitness": 0.0,
                "generation": 0,
                "max_generations": 100,
                "best_coverage": 0.0,
                "best_edge_score": 0.0,
                "theoretical_max": 0.0,
                "current_test_name": "N/A",
                "difficulty": "Unknown"  # Dodane pole difficulty
            }
            for pid in range(num_processes)
        }
        self.completed_tests = 0
        self.screen: Optional[curses.window] = None
        self.window_lock = threading.Lock()
        self.start_time = datetime.now()
        self.best_fitness_overall = float('-inf')  # Zmienione na -inf żeby obsłużyć ujemne wartości
        self.total_errors = 0

        if self.use_ui:
            self.ui_thread = threading.Thread(target=self._start_curses_ui, daemon=True)
            self.ui_thread.start()
            
    def _get_difficulty_from_path(self, path: str) -> str:
        """Wyciąga poziom trudności ze ścieżki pliku."""
        path = path.lower()
        if "easy" in path:
            return "Easy"
        elif "medium" in path:
            return "Medium"
        elif "hard" in path:
            return "Hard"
        return "Unknown"

    def _create_progress_bar(self, progress: float, width: int, fitness: float) -> str:
        filled_width = int(width * progress / 100)
        if fitness > self.best_fitness_overall:
            self.best_fitness_overall = fitness
        bar = '█' * filled_width + '░' * (width - filled_width)
        # Wyświetlamy fitness bez modyfikacji (może być ujemne)
        return f"[{bar}] {progress:5.1f}% (Best: {self.best_fitness_overall:.2f})"

    def _start_curses_ui(self) -> None:
        try:
            self.screen = curses.initscr()
            curses.start_color()
            curses.use_default_colors()
            curses.init_pair(1, curses.COLOR_RED, -1)
            curses.init_pair(2, curses.COLOR_GREEN, -1)
            curses.init_pair(3, curses.COLOR_CYAN, -1)
            curses.init_pair(4, curses.COLOR_YELLOW, -1)
            curses.init_pair(5, curses.COLOR_MAGENTA, -1)
            curses.noecho()
            curses.cbreak()
            self.screen.nodelay(True)

            self._draw_initial_screen()

            while True:
                with self.window_lock:
                    self._update_process_windows()
                    self._update_stats_window()
                    self._update_output_window()
                    self._update_footer()
                time.sleep(0.1)

        except Exception as e:
            logging.exception(f"Exception in curses UI: {e}")
        finally:
            if self.use_ui and self.screen:
                curses.nocbreak()
                self.screen.keypad(False)
                curses.echo()
                curses.endwin()

    def _draw_initial_screen(self) -> None:
        self.screen.clear()
        height, width = self.screen.getmaxyx()

        process_win_height = 4
        stats_win_height = 6
        output_win_height = height - (process_win_height * self.num_processes) - stats_win_height - 4

        process_windows = []
        for pid in range(self.num_processes):
            win = curses.newwin(process_win_height, width, pid * process_win_height, 0)
            win.attrset(curses.color_pair(3))
            win.box()
            win.addstr(0, 2, f" Process {pid} ", curses.color_pair(4) | curses.A_BOLD)
            process_windows.append(win)
            win.refresh()
        self.process_windows = process_windows

        self.stats_window = curses.newwin(
            stats_win_height,
            width,
            self.num_processes * process_win_height,
            0
        )
        self.stats_window.attrset(curses.color_pair(5))
        self.stats_window.box()
        self.stats_window.addstr(0, 2, " Global Statistics ", curses.color_pair(4) | curses.A_BOLD)
        self.stats_window.refresh()

        self.output_window = curses.newwin(
            output_win_height,
            width,
            self.num_processes * process_win_height + stats_win_height,
            0
        )
        self.output_window.attrset(curses.color_pair(2))
        self.output_window.box()
        self.output_window.addstr(0, 2, " Output Log ", curses.color_pair(4) | curses.A_BOLD)
        self.output_window.scrollok(True)
        self.output_window.refresh()

        self.footer_window = curses.newwin(2, width, height - 2, 0)
        self.footer_window.attrset(curses.color_pair(3))
        self.footer_window.box()
        self.footer_window.addstr(0, 2, " Execution Status ", curses.color_pair(4) | curses.A_BOLD)
        self.footer_window.refresh()

    def _update_process_windows(self) -> None:
        for pid, status in self.process_status.items():
            win = self.process_windows[pid]
            win.erase()
            win.box()
            difficulty_str = f" Process {pid} [{status['difficulty']}] "  # Dodane wyświetlanie difficulty
            win.addstr(0, 2, difficulty_str, curses.color_pair(4) | curses.A_BOLD)

            progress_bar = self._create_progress_bar(status['progress'], 40, status['fitness'])
            win.addstr(1, 2, progress_bar)

            gen_info = f"Gen: {status['generation']}/{status['max_generations']}"
            win.addstr(1, 55, gen_info)
            fit_info = f"Fitness: {status['fitness']:.2f}"
            win.addstr(1, 85, fit_info)

            metrics = (
                f"Cov: {status['best_coverage']:.1f} | "
                f"Edge: {status['best_edge_score']:.1f} | "
                f"MaxFit: {status['theoretical_max']:.1f}"
            )
            win.addstr(2, 2, metrics, curses.color_pair(5))

            if status['current_test_name'] != "N/A":
                win.addstr(2, 60, f"Test: {status['current_test_name']}", curses.color_pair(3))

            win.refresh()

    def _update_stats_window(self) -> None:
        self.stats_window.erase()
        self.stats_window.box()
        self.stats_window.addstr(0, 2, " Global Statistics ", curses.color_pair(4) | curses.A_BOLD)

        avg_fitness = sum(s['fitness'] for s in self.process_status.values()) / len(self.process_status)
        max_fitness = max(s['fitness'] for s in self.process_status.values())
        avg_progress = sum(s['progress'] for s in self.process_status.values()) / len(self.process_status)

        stats = [
            f"Best Fitness Overall: {self.best_fitness_overall:.2f}",
            f"Current Avg Fitness: {avg_fitness:.2f}",
            f"Current Max Fitness: {max_fitness:.2f}",
            f"Average Progress: {avg_progress:.1f}%",
            f"Error Count: {self.total_errors}"
        ]

        for i, stat in enumerate(stats, 1):
            self.stats_window.addstr(i, 2, stat, curses.color_pair(5))

        self.stats_window.refresh()

    def _update_output_window(self) -> None:
        with self.output_lock:
            for message in self.output_messages:
                if self.use_ui:
                    try:
                        if "error" in message.lower() or "exception" in message.lower():
                            self.total_errors += 1
                            self.output_window.addstr(message + "\n", curses.color_pair(1))
                        else:
                            self.output_window.addstr(message + "\n", curses.color_pair(2))
                    except curses.error:
                        self.output_window.scroll(1)
                        try:
                            self.output_window.addstr(message + "\n")
                        except curses.error:
                            pass
            self.output_messages.clear()
            if self.use_ui:
                self.output_window.refresh()

    def _update_footer(self) -> None:
        if self.use_ui:
            elapsed_time = datetime.now() - self.start_time
            elapsed_str = str(timedelta(seconds=int(elapsed_time.total_seconds())))
            progress = (self.completed_tests / self.total_tests) * 100 if self.total_tests > 0 else 0
            rate = self.completed_tests / max(elapsed_time.total_seconds() / 60, 0.001)

            self.footer_window.erase()
            self.footer_window.box()
            self.footer_window.addstr(0, 2, " Execution Status ", curses.color_pair(4) | curses.A_BOLD)

            status = (
                f"Tests: {self.completed_tests}/{self.total_tests} ({progress:.1f}%) | "
                f"Time: {elapsed_str} | "
                f"Rate: {rate:.1f} tests/min | "
                f"Active Processes: {sum(1 for s in self.process_status.values() if s['progress'] < 100)}"
            )

            self.footer_window.addstr(1, 2, status)
            self.footer_window.refresh()

    def add_output(self, message: str, error: bool = False) -> None:
        with self.output_lock:
            self.output_messages.append(message)

    def update_process(self, pid: int, progress: float, status: str, fitness: float,
                       coverage: float = 0.0, edge_score: float = 0.0, theoretical_max: float = 0.0) -> None:
        if pid not in self.process_status:
            return
        try:
            generation = 0
            max_generation = 100

            if "Gen " in status:
                part = status.split("Gen ")[1]
                gen_part = part.split()[0]
                if "/" in gen_part:
                    generation = int(gen_part.split("/")[0])
                    max_generation = int(gen_part.split("/")[1])

            self.process_status[pid].update({
                "progress": progress,
                "status": status,
                "fitness": fitness,  # Zachowujemy oryginalną wartość fitness
                "generation": generation,
                "max_generations": max_generation,
                "best_coverage": coverage,
                "best_edge_score": edge_score,
                "theoretical_max": theoretical_max
            })
        except (IndexError, ValueError) as e:
            logging.error(f"Error updating process status: {e}")
            self.process_status[pid].update({
                "progress": progress,
                "status": status,
                "fitness": fitness,
                "best_coverage": coverage,
                "best_edge_score": edge_score,
                "theoretical_max": theoretical_max
            })

    def update_test_info(self, pid: int, test_path: str) -> None:
        if pid in self.process_status:
            difficulty = self._get_difficulty_from_path(test_path)
            self.process_status[pid]["difficulty"] = difficulty
            self.process_status[pid]["current_test_name"] = test_path

    def increment_completed(self) -> None:
        self.completed_tests += 1

    def cleanup(self) -> None:
        if self.use_ui:
            try:
                self.use_ui = False
                time.sleep(0.2)
                if self.screen:
                    self.screen.keypad(False)
                    curses.nocbreak()
                    curses.echo()
                    curses.endwin()
            except Exception as e:
                logging.exception(f"Exception during cleanup: {e}")

    def wait_for_ui(self) -> None:
        if self.use_ui and hasattr(self, 'ui_thread'):
            self.ui_thread.join()
