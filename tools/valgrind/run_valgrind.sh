#!/bin/bash

# Ścieżka do wykonywalnego pliku programu w build_valgrind
EXECUTABLE="./build_valgrind/optymalizacja_kombinatoryczna"

# Funkcja do uruchamiania Valgrind
run_valgrind() {
    MODE=$1
    LOG_FILE=$2
    CMD=$3

    echo "Uruchamianie Valgrind dla trybu: $MODE"
    valgrind --tool=memcheck \
             --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             --verbose \
             --log-file=$LOG_FILE \
             $EXECUTABLE $CMD
    echo "Zakończono tryb: $MODE. Wyniki zapisane w $LOG_FILE"
    echo "---------------------------------------------"
}

# Tworzenie katalogu na logi jeśli nie istnieje
mkdir -p valgrind_logs

# Generate instance
run_valgrind "generate_instance" "valgrind_logs/valgrind_generate_instance.log" \
"generate_instance -n 400 -k 8 -dk 1 -ln 10 -lp 10 -o generated_instance.txt"

# Test instance
run_valgrind "test_instance" "valgrind_logs/valgrind_test_instance.log" \
"test_instance -i generated_instance.txt -o test_results.txt -pid 12345 -diff Medium"

# Debug mode
run_valgrind "debug" "valgrind_logs/valgrind_debug.log" \
"debug --debug"

# Tuning
run_valgrind "tuning" "valgrind_logs/valgrind_tuning.log" \
"tuning -cfg config.cfg -out tuning_results.csv"

# Tuning hybrid
run_valgrind "tuning_hybrid" "valgrind_logs/valgrind_tuning_hybrid.log" \
"tuning_hybrid -cfg config.cfg -out hybrid_tuning_results.csv"

# Combine all logs into one file
echo "Łączenie wszystkich logów..."
cat valgrind_logs/valgrind_*.log > valgrind_logs/combined_valgrind.log
echo "Wszystkie logi połączone w valgrind_logs/combined_valgrind.log" 