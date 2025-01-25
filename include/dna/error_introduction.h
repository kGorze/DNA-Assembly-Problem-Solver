#pragma once

#include "utils/logging.h"
#include "dna/dna_instance.h"
#include "utils/random.h"
#include <memory>
#include <random>
#include <mutex>
#include <algorithm>
#include <vector>
#include <string>

// Interface for error introduction strategies
class IErrorIntroductionStrategy {
public:
    virtual ~IErrorIntroductionStrategy() = default;
    virtual void introduceErrors(DNAInstance& instance) = 0;
};

// Context for error introduction
struct ErrorContext {
    int lNeg = 0;
    int lPoz = 0;
    bool repAllowed = false;
    int probablePositive = 0;
};

// Base class for error introducers
class BaseErrorIntroducer : public IErrorIntroductionStrategy {
protected:
    mutable std::mutex m_mutex;
    
    // Validate instance before introducing errors
    bool validateInstance(const DNAInstance& instance) const;
};

// Strategy for introducing negative errors (removing k-mers)
class NegativeErrorIntroducer : public BaseErrorIntroducer {
private:
    int m_lNeg;
    std::unique_ptr<Random> m_random;

public:
    explicit NegativeErrorIntroducer(int lNeg) : m_lNeg(lNeg) {
        m_random = std::make_unique<Random>();
    }
    // Add default constructor for testing
    NegativeErrorIntroducer() : NegativeErrorIntroducer(2) {}

    void introduceErrors(DNAInstance& instance) override;
};

// Strategy for introducing positive errors (adding incorrect k-mers)
class PositiveErrorIntroducer : public BaseErrorIntroducer {
private:
    int m_lPoz;
    int m_k;
    std::unique_ptr<Random> m_random;

public:
    PositiveErrorIntroducer(int lPoz, int k) : m_lPoz(lPoz), m_k(k) {
        m_random = std::make_unique<Random>();
    }
    // Add default constructor for testing
    PositiveErrorIntroducer() : PositiveErrorIntroducer(2, 5) {}

    void introduceErrors(DNAInstance& instance) override;
    std::string generateRandomKmer(int length) const;
};

// Factory for creating error introducers
class ErrorIntroducerFactory {
public:
    static std::unique_ptr<IErrorIntroductionStrategy> createNegativeErrorIntroducer(int lNeg);
    static std::unique_ptr<IErrorIntroductionStrategy> createPositiveErrorIntroducer(int lPoz, int k);
}; 