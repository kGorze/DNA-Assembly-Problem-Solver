//
// Created by konrad_guest on 08/01/2025.
//

#include "utils/utility_functions.h"


int levenshteinDistance(const std::string &s1, const std::string &s2) {
    int len1 = (int)s1.size();
    int len2 = (int)s2.size();
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1, 0));

    for(int i = 0; i <= len1; ++i) dp[i][0] = i;
    for(int j = 0; j <= len2; ++j) dp[0][j] = j;

    for(int i = 1; i <= len1; ++i) {
        for(int j = 1; j <= len2; ++j) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
            dp[i][j] = std::min({ 
                dp[i-1][j] + 1, 
                dp[i][j-1] + 1, 
                dp[i-1][j-1] + cost
            });
        }
    }
    return dp[len1][len2];
}