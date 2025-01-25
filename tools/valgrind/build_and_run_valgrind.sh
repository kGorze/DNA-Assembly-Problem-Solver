#!/bin/bash

# Kolory do wyświetlania
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Ścieżki
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." &> /dev/null && pwd)"

log() {
    echo -e "${GREEN}[BUILD+VALGRIND]${NC} $1"
}

section() {
    echo -e "\n${BLUE}=== $1 ===${NC}"
}

# Przejście do katalogu głównego projektu
cd "$PROJECT_ROOT" || exit 1

# 1. Budowanie projektu w trybie debug
section "BUDOWANIE PROJEKTU"
./tools/build.sh --debug

if [ $? -ne 0 ]; then
    echo -e "${RED}Błąd podczas budowania projektu${NC}"
    exit 1
fi

# 2. Uruchomienie Valgrinda
section "URUCHAMIANIE VALGRIND"
cd tools/valgrind || exit 1
./run_valgrind.sh

if [ $? -ne 0 ]; then
    echo -e "${RED}Błąd podczas uruchamiania Valgrind${NC}"
    exit 1
fi

section "ZAKOŃCZONO"
log "Budowanie i analiza Valgrind zakończone pomyślnie"
log "Logi Valgrind znajdują się w: tools/valgrind/valgrind_logs/" 