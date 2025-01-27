# Optymalizacja Kombinatoryczna - DNA Assembly Problem Solver

## Spis treści
- [Wprowadzenie](#wprowadzenie)
- [Struktura projektu](#struktura-projektu)
- [Algorytmy](#algorytmy)
- [Komponenty](#komponenty)
- [Konfiguracja](#konfiguracja)
- [Interfejsy](#interfejsy)
- [Szczegółowa dokumentacja](#szczegółowa-dokumentacja)

## Wprowadzenie
Program implementuje rozwiązanie problemu składania DNA (DNA Assembly Problem) przy użyciu algorytmów metaheurystycznych, w szczególności algorytmu genetycznego. Program zawiera zaawansowane mechanizmy optymalizacji i adaptacji parametrów.

## Struktura projektu
```
.
├── include/                 # Pliki nagłówkowe
│   ├── metaheuristics/     # Implementacje metaheurystyk
│   ├── interfaces/         # Interfejsy
│   ├── dna/               # Obsługa instancji DNA
│   ├── utils/             # Narzędzia pomocnicze
│   └── configuration/     # Konfiguracja
├── src/                    # Pliki źródłowe
└── test/                  # Testy jednostkowe
```

## Algorytmy

### Algorytm Genetyczny (GeneticAlgorithm)
Główna klasa implementująca algorytm genetyczny dla problemu DNA Assembly.

#### Główne komponenty:
- **Inicjalizacja**: `initializePopulation(int popSize, const DNAInstance& instance)`
  - Tworzy początkową populację rozwiązań
  - Parametry:
    - `popSize`: Rozmiar populacji
    - `instance`: Instancja problemu DNA

- **Ewaluacja**: `evaluatePopulation()`
  - Ocenia jakość każdego osobnika w populacji
  - Wykorzystuje funkcję przystosowania opartą na grafie

- **Selekcja**: Implementowana przez interfejs `ISelection`
  - Wybiera osobniki do reprodukcji
  - Dostępne strategie:
    - Selekcja turniejowa
    - Selekcja rankingowa
    - Selekcja ruletki

- **Krzyżowanie**: Implementowane przez interfejs `ICrossover`
  - Tworzy nowe rozwiązania przez kombinację istniejących
  - Dostępne operatory:
    - Krzyżowanie jednopunktowe
    - Krzyżowanie dwupunktowe
    - Krzyżowanie adaptacyjne (AdaptiveCrossover)

- **Mutacja**: Implementowana przez interfejs `IMutation`
  - Wprowadza losowe zmiany w rozwiązaniach
  - Dostępne operatory:
    - Mutacja przez zamianę
    - Mutacja przez inwersję
    - Mutacja przez przesunięcie

### Reprezentacja rozwiązań
Implementowana w klasie `PermutationRepresentation`:
- Reprezentacja permutacyjna fragmentów DNA
- Operacje na permutacjach
- Walidacja poprawności rozwiązań

## Komponenty

### DNAInstance
Klasa reprezentująca instancję problemu:
```cpp
class DNAInstance {
    std::vector<std::string> fragments;  // Fragmenty DNA
    int k;                              // Parametr k dla odległości Levenshteina
    std::string originalSequence;       // Oryginalna sekwencja (jeśli znana)
};
```

### Individual
Reprezentacja pojedynczego osobnika w populacji:
```cpp
class Individual {
    std::vector<int> genes;            // Reprezentacja rozwiązania
    double fitness;                    // Wartość funkcji przystosowania
    std::string dna;                   // Złożona sekwencja DNA
};
```

### FitnessImpl
Implementacje funkcji przystosowania:
- `OptimizedGraphBasedFitness`: Ocena bazująca na grafie
- `LevenshteinBasedFitness`: Ocena bazująca na odległości Levenshteina

## Konfiguracja

### GAConfig
Konfiguracja algorytmu genetycznego:
```cpp
struct GAConfig {
    int populationSize;                // Rozmiar populacji
    int maxGenerations;               // Maksymalna liczba generacji
    double crossoverProbability;      // Prawdopodobieństwo krzyżowania
    double mutationProbability;       // Prawdopodobieństwo mutacji
    SelectionType selectionType;      // Typ selekcji
    CrossoverType crossoverType;      // Typ krzyżowania
    MutationType mutationType;        // Typ mutacji
};
```

## Interfejsy

### IAlgorithm
```cpp
class IAlgorithm {
    virtual std::string run(const DNAInstance& instance) = 0;
    virtual ~IAlgorithm() = default;
};
```

### ISelection
```cpp
class ISelection {
    virtual std::vector<std::shared_ptr<Individual>> select(
        const std::vector<std::shared_ptr<Individual>>& population,
        int count) = 0;
};
```

### ICrossover
```cpp
class ICrossover {
    virtual std::pair<std::shared_ptr<Individual>, std::shared_ptr<Individual>>
    crossover(const std::shared_ptr<Individual>& parent1,
             const std::shared_ptr<Individual>& parent2) = 0;
};
```

### IMutation
```cpp
class IMutation {
    virtual void mutate(std::shared_ptr<Individual> individual) = 0;
};
```

## Szczegółowa dokumentacja

### Klasa GeneticAlgorithm

#### Metody publiczne:
- `GeneticAlgorithm(std::unique_ptr<IRepresentation> representation, GAConfig& config, bool debugMode = false)`
  - Konstruktor inicjalizujący algorytm genetyczny
  - Parametry:
    - `representation`: Implementacja reprezentacji rozwiązań
    - `config`: Konfiguracja algorytmu
    - `debugMode`: Tryb debugowania

- `std::string run(const DNAInstance& instance)`
  - Uruchamia algorytm genetyczny
  - Zwraca: Najlepsze znalezione rozwiązanie (sekwencję DNA)

- `double getBestFitness()`
  - Zwraca wartość przystosowania najlepszego rozwiązania

#### Metody prywatne:
- `void initializePopulation(int popSize, const DNAInstance& instance)`
  - Inicjalizuje populację początkową

- `void evaluatePopulation()`
  - Ocenia wszystkie osobniki w populacji

- `bool updateGlobalBest()`
  - Aktualizuje najlepsze znalezione rozwiązanie

- `void logGenerationStats()`
  - Zapisuje statystyki generacji

### Klasa AdaptiveCrossover

Implementuje adaptacyjny operator krzyżowania:
```cpp
class AdaptiveCrossover : public ICrossover {
    void updateCrossoverRates();      // Aktualizacja współczynników
    void adaptParameters();           // Adaptacja parametrów
    double calculateSuccess();        // Obliczenie skuteczności
};
```

### Klasa PathAnalyzer

Analizuje ścieżki w grafie reprezentującym problem:
```cpp
class PathAnalyzer {
    std::vector<PreprocessedEdge> findBestPath();
    double evaluatePath(const std::vector<int>& path);
    void preprocessEdges();
};
```

## Użycie programu

### Kompilacja
```bash
mkdir build
cd build
cmake ..
make
```

### Uruchomienie
```bash
./dna_assembly <plik_instancji> <plik_konfiguracji>
```

### Format pliku instancji
```
n k
fragment_1
fragment_2
...
fragment_n
```
gdzie:
- n: liczba fragmentów
- k: parametr k dla odległości Levenshteina
- fragment_i: i-ty fragment DNA

### Format pliku konfiguracji
```json
{
  "population_size": 100,
  "max_generations": 1000,
  "crossover_probability": 0.8,
  "mutation_probability": 0.1,
  "selection_type": "tournament",
  "crossover_type": "adaptive",
  "mutation_type": "swap"
}
```

## Wymagania systemowe
- C++17 lub nowszy
- CMake 3.10 lub nowszy
- Boost (dla testów jednostkowych)