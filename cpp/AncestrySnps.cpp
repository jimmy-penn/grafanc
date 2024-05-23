#include "AncestrySnps.h"
#include <fstream>

AncestrySnp::AncestrySnp(int id, int rsNum, int ch, int g37, int g38, char a1, char a2, float* refPops, float *vtxPops)
{
    snpId = id;
    rs = rsNum;
    chr = ch;
    posG37 = g37;
    posG38 = g38;
    ref = a1;
    alt = a2;

    for (int i = 0; i < numRefPops; i++) refPopAfs[i] = refPops[i];
    for (int i = 0; i < numVtxPops; i++) vtxPopAfs[i] = vtxPops[i];
}

void AncestrySnp::SetRefSubPopAfs (float* popAfs)
{
    // By default, AFs for normalization is the same as reference AFs
    for (int i = 0; i < numSubPops; i++) {
        refSubPopAfs[i] = popAfs[i];
        nomSubPopAfs[i] = popAfs[i];
    }
}

void AncestrySnp::SetNomSubPopAf (int popId, float popAf)
{
    // When needed, AF for normalization can be set to different value then ref AF
    nomSubPopAfs[popId] = popAf;
}

AncestrySnps::AncestrySnps()
{
  
  refPopNames[0] = "African";
  refPopNames[1] = "European";
  refPopNames[2] = "Asian";
  refPopNames[3] = "Mexican";
  refPopNames[4] = "Indian-Pakistani";
  
  string subPopNames[numSubPops];

  subPopNames[0] = "chn"; // Chinease, especially Northern Chinese
  subPopNames[1] = "jpn"; // Japanese
  subPopNames[2] = "sea"; // Southeast Asian: Thai, Vietnamese, ...
  subPopNames[3] = "pac"; // Pacific Islander
  subPopNames[4] = "ceu"; // English, Irish, ...
  subPopNames[5] = "seu"; // South European: Italian, ...
  subPopNames[6] = "fin"; // Finnish
  subPopNames[7] = "lt1"; // Lating American 1
  
  for (int popId = 0; popId < numSubPops; popId++) {
      refSubPopNames[popId] = "ref_" + subPopNames[popId];
      nomSubPopNames[popId] = "nom_" + subPopNames[popId];
  }

  snps = {};
}


AncestrySnps::~AncestrySnps()
{
    snps.clear();
    rsToAncSnpId.clear();
    pos37ToAncSnpId.clear();
    pos38ToAncSnpId.clear();
}

int AncestrySnps::ReadAncestrySnpsFromFile(string ancSnpFile)
{
    ASSERT(FileExists(ancSnpFile.c_str()), "File " << ancSnpFile << " does not exist!\n");

    double popExpPfSums[numRefPops];
    double popExpPaSums[numRefPops];
    double popExpPeSums[numRefPops];

    for (int popId = 0; popId < numRefPops; popId++) {
        popExpPeSums[popId] = 0;
        popExpPfSums[popId] = 0;
        popExpPaSums[popId] = 0;
    }

    int lineLen = 300;
    char snpLine[lineLen];

    FILE *ifp = fopen(ancSnpFile.c_str(), "r");
    ASSERT(ifp, "Couldn't open file " << ancSnpFile << ".\n");

    int lineNo = 0;
    bool fileIsValid = true;

    int snpId = 0;
    int rsNum, chr, g37, g38;
    char a1, a2;
    float rfEur, rfAfa, rfAsn, rfLat, rfSas, vtEur, vtAfr, vtEas;

    while (fgets(snpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        if (lineNo == 0) {
            if (snpLine[0] != 'c' || snpLine[1] != 'h' || snpLine[2] != 'r') {
                fileIsValid = false;
            }
        }
        else {
            sscanf(snpLine, "%d %d %d %d %c %c %f %f %f %f %f %f %f %f",
            &chr, &g37, &g38, &rsNum, &a1, &a2, &rfEur, &rfAfa, &rfAsn, &rfLat, &rfSas, &vtEur, &vtAfr, &vtEas);

            float* refPopAfs = new float[numRefPops];
            float* vtxPopAfs = new float[numVtxPops];

            refPopAfs[0] = rfEur;
            refPopAfs[1] = rfAfa;
            refPopAfs[2] = rfAsn;
            refPopAfs[3] = rfLat;
            refPopAfs[4] = rfSas;

            vtxPopAfs[0] = vtEur;
            vtxPopAfs[1] = vtAfr;
            vtxPopAfs[2] = vtEas;

            AncestrySnp ancSnp(snpId, rsNum, chr, g37, g38, a1, a2, refPopAfs, vtxPopAfs);

            snps.push_back(ancSnp);

            long int chrPos37 = (long)chr * 1000000000 + g37;
            long int chrPos38 = (long)chr * 1000000000 + g38;

            rsToAncSnpId[rsNum] = snpId;
            pos37ToAncSnpId[chrPos37] = snpId;
            pos38ToAncSnpId[chrPos38] = snpId;

            double pev = refPopAfs[0];
            double pfv = refPopAfs[1];
            double pav = refPopAfs[2];

            double qev = 1 - pev;
            double qfv = 1 - pfv;
            double qav = 1 - pav;

            for (int vtxId = 0; vtxId < 3; vtxId++) {
                double pv  = vtxPopAfs[vtxId];
                double qv  = 1 - pv;

                double aaPev = log(pev) * 2;
                double bbPev = log(qev) * 2;
                double abPev = log(pev) + log(qev) + log(2);

                double aaPfv = log(pfv) * 2;
                double bbPfv = log(qfv) * 2;
                double abPfv = log(pfv) + log(qfv) + log(2);

                double aaPav = log(pav) * 2;
                double bbPav = log(qav) * 2;
                double abPav = log(pav) + log(qav) + log(2);

                double eGd = aaPev * pv * pv + bbPev * qv * qv + abPev * 2 * pv * qv;
                double fGd = aaPfv * pv * pv + bbPfv * qv * qv + abPfv * 2 * pv * qv;
                double aGd = aaPav * pv * pv + bbPav * qv * qv + abPav * 2 * pv * qv;

                vtxExpGenoDists[vtxId][0][snpId] = eGd;
                vtxExpGenoDists[vtxId][1][snpId] = fGd;
                vtxExpGenoDists[vtxId][2][snpId] = aGd;

                popExpPeSums[vtxId] += eGd;
                popExpPfSums[vtxId] += fGd;
                popExpPaSums[vtxId] += aGd;
            }

            delete refPopAfs;
            delete vtxPopAfs;

            snpId++;
        }

        lineNo++;
    }
    fclose(ifp);

    ASSERT(snpId == numAllAncSnps, "snpId = " << numAllAncSnps << ".\n");

    for (int vtxId = 0; vtxId < 3; vtxId++) {
        vtxPopExpGds[vtxId].e = -1 * popExpPeSums[vtxId]/snpId;
        vtxPopExpGds[vtxId].f = -1 * popExpPfSums[vtxId]/snpId;
        vtxPopExpGds[vtxId].a = -1 * popExpPaSums[vtxId]/snpId;
    }

    cout << "Read " << snpId << " ancestry SNPs from file " << ancSnpFile << "\n\n";
    if (0) {
        cout << "Expected vertex genetic distances\n";
        for (int vtxId = 0; vtxId < 3; vtxId++) {
            cout << "\tVertex " << vtxId << "\n";
            cout << "\t\tEUR: " << vtxPopExpGds[vtxId].e << "\n";
            cout << "\t\tAFR: " << vtxPopExpGds[vtxId].f << "\n";
            cout << "\t\tEAS: " << vtxPopExpGds[vtxId].a << "\n";
        }
    }
    
    return snpId;
}

int AncestrySnps::ReadRefSubPopSnpsFromFile(string refPopFile)
{
    ASSERT(FileExists(refPopFile.c_str()), "File " << refPopFile << " does not exist!\n");

    cout << "Reading AFs from file " << refPopFile << "\n";  
    
    int snpId = 0;
    
    FILE* fp = fopen(refPopFile.c_str(), "r");
    if (fp == NULL) exit(EXIT_FAILURE);
    
    char* line = NULL;
    size_t len = 0;
    const string delim = "\t";

    // Make sure the file has correct columns
    bool isRightFile = true;
    if ((getline(&line, &len, fp)) != -1) {
        const string header(line);
        vector<string> vars = SplitString(header, delim);
        int numVars = vars.size();

        if (numVars <= numSubPops) {
            isRightFile = false;
        }
        else if (vars[0] != "rs") {
            isRightFile = false;
        }
        else {
            for (int i = 0; i < numSubPops; i++) {
                string var = TrimString(vars[i+1]);

                if (var != refSubPopNames[i]) {
                    isRightFile = false;
                    break;
                }
            }
        }
    }

    string expCols = "rs";
    for (int i = 0; i < numSubPops; i++) {
        expCols += ", " + refSubPopNames[i];
    }
    
    ASSERT(isRightFile, "File " << refPopFile << " does not have expected columns (" << expCols << ")!\n");
    
    while ((getline(&line, &len, fp)) != -1) {
        ASSERT(snpId < numAllAncSnps, "File " << refPopFile << " has too many rows!\n");
      
        vector<string> vals = SplitString(string(line), delim);
        int rs = stoi(vals[0]); 
        ASSERT(rs == snps[snpId].rs, "ERROR in " << refPopFile << ": Row #" << snpId << " has rs" << rs << " not " << snps[snpId].rs << "!\n");
        
        float popAfs[numSubPops];
        for (int i = 0; i < numSubPops; i++) popAfs[i] = stof(vals[i+1]);
        snps[snpId].SetRefSubPopAfs(popAfs);
        
        snpId++;
    }
    fclose(fp);
    
    if (line) free(line);

    return snpId;
}

int AncestrySnps::ReadNomSubPopSnpsFromFile(string nomPopFile)
{
    ASSERT(FileExists(nomPopFile.c_str()), "File " << nomPopFile << " does not exist!\n");
    
    cout << "Reading AFs for normalization from file " << nomPopFile << "\n";  
    
    int snpId = 0;
    
    FILE* fp = fopen(nomPopFile.c_str(), "r");
    if (fp == NULL) exit(EXIT_FAILURE);
    
    char* line = NULL;
    size_t len = 0;
    const string delim = "\t";
    
    // Find columns corresponding to each ref population
    int nomPopColIds[numSubPops];
    for (int i = 0; i < numSubPops; i++) nomPopColIds[i] = 0;

    bool isRightFile = true;
    vector<string> vars;
    int numVars = 0;
    
    if ((getline(&line, &len, fp)) != -1) {
        const string header(line);
        vars = SplitString(header, delim);
        numVars = vars.size();
        
        if (numVars < 2) {
            isRightFile = false;
        }
        else if (vars[0] != "rs") {
            isRightFile = false;
        }
    }

    ASSERT(isRightFile, "File " << nomPopFile << " does not have expected columns (rs, ...)!\n");
    
    int numPopsInFile = 0;
    
    for (int i = 0; i < numSubPops; i++) {
        string expVar = nomSubPopNames[i];      
        for (int col = 1; col < numVars; col++) {
            string var = TrimString(vars[col]);

            if (var == expVar) {
                nomPopColIds[i] = col;
                numPopsInFile++;
                break;
            }
        }
    }        

    if (numPopsInFile < 1) {
        cout << "ERROR: no reference populations are found in " << nomPopFile << "\n";
    }
    else {
        while ((getline(&line, &len, fp)) != -1) {
            ASSERT(snpId < numAllAncSnps, "File " << nomPopFile << " has too many rows!\n");
            
            vector<string> vals = SplitString(string(line), delim);
            int rs = stoi(vals[0]); 
            ASSERT(rs == snps[snpId].rs, "ERROR in " << nomPopFile << ": Row #" << snpId << " has rs" << rs << " not " << snps[snpId].rs << "!\n");

            for (int i = 0; i < numSubPops; i++) {
                int col = nomPopColIds[i];

                if (col > 0) {
                    float af = stof(vals[col]);
                    snps[snpId].SetNomSubPopAf(i, af);
                }
            }
            
            snpId++;
        }
    }
    
    fclose(fp);
    
    if (line) free(line);

    return snpId;
}

int AncestrySnps::FindSnpIdGivenRs(int rsNum)
{
    int snpId = -1;

    if (rsToAncSnpId.find(rsNum) != rsToAncSnpId.end()) {
        snpId = rsToAncSnpId[rsNum];
    }

    return snpId;
}

int AncestrySnps::FindSnpIdGivenChrPos(int chr, int pos, int build)
{
    int snpId = -1;

    long int chrPos = long(chr) * 1000000000 + pos;

    if (build == 37) {
        if (pos37ToAncSnpId.find(chrPos) != pos37ToAncSnpId.end()) {
            snpId = pos37ToAncSnpId[chrPos];
        }
    }
    else if (build == 38) {
        if (pos38ToAncSnpId.find(chrPos) != pos38ToAncSnpId.end()) {
            snpId = pos38ToAncSnpId[chrPos];
        }
    }

    return snpId;
}

AncestrySnp AncestrySnps::GetAncestrySnp(int snpId)
{
    return snps[snpId];
}

void AncestrySnps::ShowAncestrySnps()
{
    int numAncSnps = snps.size();

    cout << "Total " << numAncSnps << " Ancestry SNPs.\n";
    bool debug = 0;

    if (debug) {
        for (int i = 0; i < 20; i++) {
            int snpId = i * 5000;
            AncestrySnp snp = snps[snpId];
            cout << "SNP " << snpId << " rs " << snp.rs << " chr " << snp.chr << " pos " << snp.posG37
            << " ref " << snp.ref << " alt " << snp.alt;
            for (int j = 0; j < numRefPops; j++) cout << " Ref " << j << " = " << snp.refPopAfs[j];
            for (int j = 0; j < numVtxPops; j++) cout << " Vtx " << j << " = " << snp.vtxPopAfs[j];
            cout << "\n";
        }
    }

    cout << "Positions (x, y, z coordinates) of the three vertices when all SNPs have genotypes:\n";
    printf("\tE:  %5.4f  %5.4f  %5.4f\n", vtxPopExpGds[0].e, vtxPopExpGds[0].f, vtxPopExpGds[0].a);
    printf("\tF:  %5.4f  %5.4f  %5.4f\n", vtxPopExpGds[1].e, vtxPopExpGds[1].f, vtxPopExpGds[1].a);
    printf("\tA:  %5.4f  %5.4f  %5.4f\n", vtxPopExpGds[2].e, vtxPopExpGds[2].f, vtxPopExpGds[2].a);

    cout << "\nExpected genetic distances\n";
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << i << "-" << j << ": ";
            for (int k = 0; k < 5; k++) {
                cout << vtxExpGenoDists[i][j][k] << " ";
            }
            cout << "\n";
        }
    }
}

