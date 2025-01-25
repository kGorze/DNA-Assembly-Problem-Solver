# Optymalizacja Kombinatoryczna - Projekt

## Struktura projektu

```
.
├── src/                    # Kod źródłowy C++
├── include/               # Pliki nagłówkowe C++
├── scripts/               # Skrypty pomocnicze (Python)
│   ├── sbh_testing_framework.py
│   └── display_manager.py
├── tools/                # Narzędzia deweloperskie
│   └── valgrind/        # Konfiguracja i skrypty Valgrind
├── config/               # Pliki konfiguracyjne
├── tests/                # Testy i ich wyniki
│   ├── configs/         # Konfiguracje testów
│   └── outputs/         # Wyniki testów
├── build_outputs/        # Katalog na pliki wyjściowe kompilacji
├── google_tests/         # Testy jednostkowe
└── CMakeLists.txt       # Główny plik konfiguracyjny CMake
```

## Opis katalogów

- `src/` i `include/` - Główny kod źródłowy projektu
- `scripts/` - Skrypty pomocnicze do testowania i wizualizacji
- `tools/valgrind/` - Narzędzia do analizy pamięci z Valgrind
- `config/` - Pliki konfiguracyjne dla różnych trybów działania
- `tests/` - Wyniki testów i ich konfiguracje
- `build_outputs/` - Miejsce na pliki wyjściowe kompilacji
- `google_tests/` - Testy jednostkowe

## Kompilacja i uruchamianie

### Standardowa kompilacja
```bash
mkdir -p build && cd build
cmake ..
make
```

### Kompilacja z Valgrind
```bash
./tools/valgrind/configure_valgrind_build.ps1
```

### Uruchamianie testów z Valgrind
```bash
cd tools/valgrind
./run_valgrind.sh
```

## Wymagania
- CMake 3.10+
- Kompilator C++ z obsługą C++17
- Python 3.8+ (dla skryptów pomocniczych)
- WSL z zainstalowanym Valgrind (dla analizy pamięci)