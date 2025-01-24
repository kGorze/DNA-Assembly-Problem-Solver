class NoImprovementStopping : public IStopping {
public:
    explicit NoImprovementStopping(int maxGenerationsWithoutImprovement);
    
    bool stop(const std::vector<std::shared_ptr<std::vector<int>>>& population,
              const DNAInstance& instance,
              int generation,
              double bestFitness) const override;

    void reset() override;

private:
    int m_maxGenerationsWithoutImprovement;
    mutable double m_lastBestFitness;
    mutable int m_generationsWithoutImprovement;
}; 