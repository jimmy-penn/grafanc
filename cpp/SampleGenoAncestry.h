#ifndef SAMPLE_GENO_ANCESTRY_H
#define SAMPLE_GENO_ANCESTRY_H

#include <fstream>
#include "Util.h"
#include "AncestrySnps.h"
#include "FamFileSamples.h"
#include "SampleGenoDist.h"

static const int numSubPopScores = 17;

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
    // Values of the 8 bits, for coding 4 genotypes in one byte. The first geno is stored in the first two bits
    const unsigned int charGenoBitVals[8] = {128, 64, 32, 16, 8, 4, 2, 1};    
  
    vector<GenoSample> samples;
    vector<int> *ancSnpIds;
    vector<unsigned char*> *ancSnpCodedGenos; // Use char, instead of int, to save space
    
    //                                            0         1          2            3         4          5          6        7         8   
    const string ancSnpFileCols[numSubPops] = {"UKBBEUR", "UKBBAFR", "UKBBEAS", "Nigeria", "Ghana", "Zimbabwe", "Uganda", "Kenya", "Barbados", 
                                     
    //        9       10            11            12        13       14        15                                           
          "China", "HongKong", "Philippines", "Thailand", "Japan", "Nepal", "VietNam",
                                     
    //       16           17           18         19        20      21     22      23    24      25      26         27        28                                 
         "Pakistan", "Bangladesh", "SriLanka", "India2", "Kenya2", "CEU", "GBR", "FIN", "TSI", "PUR", "MXL_CLM", "UKBSAS", "UKBMEX"};

    const string popScoreNames[numSubPopScores] = {"EA1", "EA2", "EA3", "EA4", "EA5", "EA6", "AF1", "AF2", "AF3", "EU1", "EU2", "SA1", "SA2", "SA3", "SA4", "GL1", "GL2"};
    
    // Pair of ref pops to use for each score   
    //                                         EA1 EA2 EA3 EA4 EA5 EA6 AF1 AF2 AF3 EU1 EU2 SA1 SA2 SA3 SA4 GL1 GL2
    const int scorePopIdx1[numSubPopScores] = { 9,  9,  9,  9, 11, 11,  3,  3,  6, 21, 21, 18, 19, 18, 18, 28, 24};
    const int scorePopIdx2[numSubPopScores] = {11, 13, 12, 14, 12, 13,  5,  4,  8, 24, 23, 17, 16, 19, 16, 27, 25};
    
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
    void SetSnpGenoData(vector<int>*, vector<unsigned char*>*);
    void SetNumThreads(int);
    void InitPopPvalues();
    void CalculateSubPopGd0Values();

    int GetNumSamples() { return numSamples; };
    int GetNumAncSamples() { return numAncSmps; };
    bool HasEnoughAncestrySnps(int numSnps) { return numSnps >= minAncSnps; }

    void ShowSummary();
};

#endif
