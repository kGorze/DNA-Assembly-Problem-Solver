std::vector<std::shared_ptr<Individual>> selected;
selected.reserve(count);

LOG_DEBUG("Tournament selection - population size: " + std::to_string(population.size()) + 
          ", tournament size: " + std::to_string(m_tournamentSize) + 
          ", selecting " + std::to_string(count) + " individuals");

for (size_t i = 0; i < count; ++i) {
    // Select tournament participants randomly
    std::vector<std::shared_ptr<Individual>> tournament;
    tournament.reserve(m_tournamentSize);
    
    for (size_t j = 0; j < m_tournamentSize; ++j) {
        size_t idx = m_random.nextInt(0, population.size() - 1);
        tournament.push_back(population[idx]);
    }
    
    // Find the best individual in tournament
    auto best = tournament[0];
    for (size_t j = 1; j < tournament.size(); ++j) {
        if (tournament[j]->getFitness() > best->getFitness()) {
            best = tournament[j];
        }
    }
    
    LOG_DEBUG("Tournament winner " + std::to_string(i+1) + " - genes: " + 
              std::to_string(best->getGenes().size()) + ", fitness: " + 
              std::to_string(best->getFitness()) + ", valid: " + 
              std::to_string(best->isValid()));
              
    selected.push_back(best);
}

return selected; 