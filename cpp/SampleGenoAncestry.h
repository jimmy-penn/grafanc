#ifndef SAMPLE_GENO_ANCESTRY_H
#define SAMPLE_GENO_ANCESTRY_H

#include <fstream>
#include "Util.h"
#include "AncestrySnps.h"
#include "FamFileSamples.h"
#include "SampleGenoDist.h"

static const int numSubPopScores = 2;

class GenoSample
{
public:
    string name;
    string father;
    string mother;
    int sex;

    // Ancestry results calculated from genotypes
    int numAncSnps;
    bool ancIsSet;
    float gd1, gd2, gd3;
    float ePct, fPct, aPct;   // Ancestry (EUR, AFR, EAS) components of the sample
    float subPopScores[numSubPopScores];

public:
    GenoSample(string);
    void SetAncestryScores(int, float, float, float, float, float, float, bool);
    void SetSubPopGdScores(float*);
};


class SampleGenoAncestry
{
private:
    int numSamples;
    int numAncSmps;

    int minAncSnps;
    int totAncSnps;
    int numAncSnps;
    int numThreads;                // Number of threads for parallel computing

    AncestrySnps *ancSnps;
    SampleGenoDist *vtxExpGd0;    // Genetic distances from 3 vertices to ref populations when all SNPs have genotypes

public:
    vector<GenoSample> samples;
    vector<int> *ancSnpIds;
    vector<char*> *ancSnpCodedGenos; // Use char, instead of int, to save space
    const string popScoreNames[numSubPopScores] = {"EA1", "EA2"};
    
    // Pair of ref pops to use for each score    
    // Ref pops are: chn, jpn, sea, pac, ceu, seu, fin
    const int scorePopIdx1[numSubPopScores] = {0, 0};
    const int scorePopIdx2[numSubPopScores] = {2, 1};
    
    // Expected GD scores of subpopulations when all SNPs are included
    float subPopGd0P1[numSubPopScores];
    float subPopGd0P2[numSubPopScores];
    
    // Normalizing GD scores to subPopGd0P1 = -1 and subPopGd0P2 = 1
    const float subPopGdNormP1 = -1.0;
    const float subPopGdNormP2 = 1.0;
    
    SampleGenoAncestry(AncestrySnps*, int=100);
    ~SampleGenoAncestry();

    void SetGenoSamples(const vector<string>&);
    void SetGenoSamples(const vector<FamSample>&);
    int SaveAncestryResults(string);
    void SetAncestryPvalues(int);
    void SetSnpGenoData(vector<int>*, vector<char*>*);
    void SetNumThreads(int);
    void InitPopPvalues();
    void CalculateSubPopGd0Values();

    int GetNumSamples() { return numSamples; };
    int GetNumAncSamples() { return numAncSmps; };
    bool HasEnoughAncestrySnps(int numSnps) { return numSnps >= minAncSnps; }

    void ShowSummary();
};

#endif
