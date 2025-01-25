#!/bin/bash

# Ścieżki do plików
EXECUTABLE="../../build_valgrind/optymalizacja_kombinatoryczna"
CONFIG_DIR="../../config"
OUTPUT_DIR="../../tests/outputs"
LOGS_DIR="./valgrind_logs"

# Maksymalny czas wykonania dla każdego testu (w sekundach)
MAX_TIME=60

# Funkcja do uruchamiania Valgrind z timeoutem
run_valgrind() {
    MODE=$1
    LOG_FILE=$2
    CMD=$3

    echo "Uruchamianie Valgrind dla trybu: $MODE"
    timeout $MAX_TIME valgrind --tool=memcheck \
             --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             --verbose \
             --log-file=$LOG_FILE \
             $EXECUTABLE $CMD

    if [ $? -eq 124 ]; then
        echo "Test przerwany po $MAX_TIME sekundach"
        echo "Test przerwany po $MAX_TIME sekundach (timeout)" >> $LOG_FILE
    else
        echo "Zakończono tryb: $MODE"
    fi
    echo "Wyniki zapisane w $LOG_FILE"
    echo "---------------------------------------------"
}

# Tworzenie katalogów na logi i wyniki
mkdir -p $LOGS_DIR
mkdir -p $OUTPUT_DIR

# Generate instance (mała instancja dla szybkich testów)
run_valgrind "generate_instance" "$LOGS_DIR/valgrind_generate_instance.log" \
"generate_instance -n 100 -k 6 -dk 1 -ln 8 -lp 8 -o $OUTPUT_DIR/generated_instance.txt"

# Test instance (krótki test)
run_valgrind "test_instance" "$LOGS_DIR/valgrind_test_instance.log" \
"test_instance -i $OUTPUT_DIR/generated_instance.txt -o $OUTPUT_DIR/test_results.txt -pid 12345 -diff Easy"

# Debug mode (podstawowe testy)
run_valgrind "debug" "$LOGS_DIR/valgrind_debug.log" \
"debug --quick-test"

# Tuning (ograniczony czas i mniejsza liczba iteracji)
run_valgrind "tuning" "$LOGS_DIR/valgrind_tuning.log" \
"tuning -cfg $CONFIG_DIR/config.cfg -out $OUTPUT_DIR/tuning_results.csv -max-time 45 -pop-size 10 -generations 5"

# Tuning hybrid (ograniczony czas i mniejsza liczba iteracji)
run_valgrind "tuning_hybrid" "$LOGS_DIR/valgrind_tuning_hybrid.log" \
"tuning_hybrid -cfg $CONFIG_DIR/config.cfg -out $OUTPUT_DIR/hybrid_tuning_results.csv -max-time 45 -pop-size 10 -generations 5"

# Łączenie wszystkich logów
echo "Łączenie wszystkich logów..."
cat $LOGS_DIR/valgrind_*.log > $LOGS_DIR/combined_valgrind.log
echo "Wszystkie logi połączone w $LOGS_DIR/combined_valgrind.log"

# Podsumowanie
echo "Podsumowanie wykonanych testów:"
echo "1. Generate Instance: mała instancja (n=100, k=6)"
echo "2. Test Instance: tryb Easy"
echo "3. Debug: szybki test"
echo "4. Tuning: 45s, populacja 10, 5 generacji"
echo "5. Tuning Hybrid: 45s, populacja 10, 5 generacji" 