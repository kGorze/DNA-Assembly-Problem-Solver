#pragma once

#include <memory>
#include <vector>
#include "dna_instance.h"
#include "individual.h"

class IRepresentation {
public:
    virtual ~IRepresentation() = default;
    
    virtual bool initializeIndividual(Individual& individual, 
                                    const DNAInstance& instance) = 0;
    
    virtual bool isValid(const std::shared_ptr<Individual>& individual,
                        const DNAInstance& instance) const = 0;
                        
    virtual std::string toString(const std::shared_ptr<Individual>& individual,
                               const DNAInstance& instance) const = 0;
};

class PermutationRepresentation : public IRepresentation {
public:
    bool initializeIndividual(Individual& individual, 
                            const DNAInstance& instance) override;
    
    bool isValid(const std::shared_ptr<Individual>& individual,
                const DNAInstance& instance) const override;
                
    std::string toString(const std::shared_ptr<Individual>& individual,
                        const DNAInstance& instance) const override;
                        
private:
    bool validateGenes(const std::vector<int>& genes,
                      const DNAInstance& instance) const;
}; 