#!/bin/bash

# Kolory do wyświetlania
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Ścieżki katalogów
BUILD_DIR="cmake-build-current" # Tymczasowy katalog do budowania
RUNTIME_DIR="build"            # Katalog z plikami wykonywalnymi i zasobami
BUILD_TYPE="Debug"
JOBS=$(nproc)

# Funkcja do wyświetlania komunikatów
log() {
    echo -e "${GREEN}[BUILD]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Funkcja do bezpiecznego kopiowania plików
safe_copy() {
    local src="$1"
    local dest="$2"
    if [ -f "$src" ]; then
        if [ -f "$dest" ]; then
            warning "Plik $dest już istnieje, tworzenie kopii zapasowej..."
            cp "$dest" "${dest}.bak"
        fi
        cp "$src" "$dest"
        log "Skopiowano: $src -> $dest"
    fi
}

# Sprawdzenie argumentów
while [[ $# -gt 0 ]]; do
    case $1 in
        --release)
            BUILD_TYPE="Release"
            BUILD_DIR="cmake-build-release"
            shift
            ;;
        --debug)
            BUILD_TYPE="Debug"
            BUILD_DIR="cmake-build-debug"
            shift
            ;;
        --clean)
            rm -rf "$BUILD_DIR"
            log "Wyczyszczono katalog $BUILD_DIR"
            shift
            ;;
        --clean-all)
            rm -rf cmake-build-*
            log "Wyczyszczono wszystkie katalogi build"
            shift
            ;;
        --jobs=*)
            JOBS="${1#*=}"
            shift
            ;;
        *)
            error "Nieznany argument: $1"
            exit 1
            ;;
    esac
done

# Tworzenie katalogów jeśli nie istnieją
mkdir -p "$BUILD_DIR"
mkdir -p "$RUNTIME_DIR"

# Przejście do katalogu build
cd "$BUILD_DIR" || exit 1

# Konfiguracja CMake jeśli nie istnieje plik cache
if [ ! -f "CMakeCache.txt" ]; then
    log "Konfiguracja CMake (BUILD_TYPE=$BUILD_TYPE)..."
    cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ..
    if [ $? -ne 0 ]; then
        error "Błąd konfiguracji CMake"
        exit 1
    fi
fi

# Sprawdzenie czy są zmiany w plikach CMake
if [ -f "CMakeCache.txt" ] && [ -n "$(find .. -name 'CMakeLists.txt' -newer CMakeCache.txt)" ]; then
    warning "Wykryto zmiany w plikach CMake, rekonfiguracja..."
    cmake ..
fi

# Budowanie projektu
log "Budowanie projektu używając $JOBS wątków..."
cmake --build . --config $BUILD_TYPE -j $JOBS

if [ $? -eq 0 ]; then
    log "Budowanie zakończone sukcesem"
    
    # Kopiowanie pliku wykonywalnego do katalogu runtime
    if [ -f "optymalizacja_kombinatoryczna" ]; then
        safe_copy "optymalizacja_kombinatoryczna" "../$RUNTIME_DIR/optymalizacja_kombinatoryczna"
        log "Program dostępny w: $RUNTIME_DIR/optymalizacja_kombinatoryczna"
    fi
else
    error "Błąd podczas budowania"
    exit 1
fi

# Powrót do katalogu głównego
cd .. 