#pragma once

#include "utils/logging.h"
#include "dna/dna_instance.h"
#include <memory>
#include <random>
#include <mutex>
#include <algorithm>

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
    mutable std::mt19937 m_rng;

public:
    explicit NegativeErrorIntroducer(int lNeg);
    void introduceErrors(DNAInstance& instance) override;
};

// Strategy for introducing positive errors (adding incorrect k-mers)
class PositiveErrorIntroducer : public BaseErrorIntroducer {
private:
    int m_lPoz;
    int m_k;
    mutable std::mt19937 m_rng;

public:
    PositiveErrorIntroducer(int lPoz, int k);
    void introduceErrors(DNAInstance& instance) override;
    std::string generateRandomKmer(int length) const;
};

// Factory for creating error introducers
class ErrorIntroducerFactory {
public:
    static std::unique_ptr<IErrorIntroductionStrategy> createNegativeErrorIntroducer(int lNeg);
    static std::unique_ptr<IErrorIntroductionStrategy> createPositiveErrorIntroducer(int lPoz, int k);
}; 