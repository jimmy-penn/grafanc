#include "VcfSampleAncestrySnpGeno.h"

VcfSampleAncestrySnpGeno::VcfSampleAncestrySnpGeno(string file, AncestrySnps *aSnps)
{
    ancSnps = aSnps;
    vcfFile = file;

    vcfSamples = {};
    vcfAncSnpGtVals = {};
    vcfAncSnpChrs = {};
    vcfAncSnpPoss = {};
    vcfAncSnpSnps = {};
    vcfAncSnpRefs = {};
    vcfAncSnpAlts = {};
    vcfRsIdAncSnpIds = {};
    vcfGb37AncSnpIds = {};
    vcfGb38AncSnpIds = {};
    vcfAncSnpCodedGenos = {};
    vcfAncSnpIds = {};

    totAncSnps = ancSnps->GetNumAncestrySnps();
    numSamples = 0;
    totVcfSnps = 0;
    putativeAncSnps = 0;
    numRsIdAncSnps = 0;
    numGb37AncSnps = 0;
    numGb38AncSnps = 0;
    numVcfAncSnps = 0;
}

VcfSampleAncestrySnpGeno::~VcfSampleAncestrySnpGeno()
{
    vcfSamples.clear();
    vcfAncSnpChrs.clear();
    vcfAncSnpPoss.clear();
    vcfAncSnpSnps.clear();
    vcfAncSnpRefs.clear();
    vcfAncSnpAlts.clear();
    vcfRsIdAncSnpIds.clear();
    vcfGb37AncSnpIds.clear();
    vcfGb38AncSnpIds.clear();
    vcfAncSnpIds.clear();

    DeleteAncSnpCodedGenos();
    DeleteAncSnpGtValues();
}

void VcfSampleAncestrySnpGeno::DeleteAncSnpCodedGenos()
{
    for (int i = 0; i < vcfAncSnpCodedGenos.size(); i++) {
        if (vcfAncSnpCodedGenos[i]) delete vcfAncSnpCodedGenos[i];
    }
    vcfAncSnpCodedGenos.clear();
}

void VcfSampleAncestrySnpGeno::DeleteAncSnpGtValues()
{
    for (int i = 0; i < vcfAncSnpGtVals.size(); i++) {
        if (!vcfAncSnpGtVals[i].empty()) vcfAncSnpGtVals[i].clear();
    }
    vcfAncSnpGtVals.clear();
}

bool VcfSampleAncestrySnpGeno::ReadDataFromFile()
{
    cout << "Reading data from file " << vcfFile << "\n";
    gzFile file = gzopen (vcfFile.c_str(), "r");

    unsigned int charBaseInts[4] = {64, 16, 4, 1};
    
    int lineNo = 0;
    int numVcfSnps = 0;
    char colValue[WORDLEN];

    vector<int> ancSnpIds;
    vector<string> chromosomes;
    vector<string> positions;
    vector<string> ids;
    vector<string> refs;
    vector<string> alts;

    colValue[0] = 0;
    int valPos = 0;
    bool fileDone = false;
    bool hasHeadRow = false;
    int buffNo = 0;
    int maxLineLen = 0;
    int vcfColNo = 0;
    int numSnpLines = 0;

    string gtyStr, chrStr, posStr, snpStr, refStr, altStr;

    vector<string> snpGts;

    while (!fileDone) {
        int err;
        int bytesRead;
        char buffer[BUFFERLEN];
        bytesRead = gzread(file, buffer, BUFFERLEN);
        buffer[bytesRead] = '\0';

        if (bytesRead < BUFFERLEN - 1) {
            if (gzeof(file)) {
                fileDone = true;
            }
            else {
                const char * errorString;
                errorString = gzerror(file, & err);
                if (err) {
                    fprintf(stderr, "Error: %s.\n", errorString);
                    return false;
                }
            }
        }

        // Read vcf file line by line, finding sample set from #CHROM row, then reading SNP genotype from the rest rows.
        // Different from the PLINK set, from which ancestry SNPs and their alleles can be determined by checking the .bim file,
        // vcf file stores all information in one file. We have to read and save the genotypes of potential ancestry SNPs first,
        // and determine which genome build to use and flip the alleles later, in order to avoid opening the large file again.
        int buffPos = 0;
        while (buffPos < bytesRead) {
            bool isNewLine = false;
            if (buffer[buffPos] == '\t' || buffer[buffPos] == '\n') {
                if (buffer[buffPos] == '\n') isNewLine = true;

                colValue[valPos] = 0;

                if      (vcfColNo == 0) {
                    chrStr = string(colValue);
                    chromosomes.push_back(chrStr);
                }
                else if (vcfColNo == 1) {
                    posStr = string(colValue);
                    positions.push_back(posStr);
                }
                else if (vcfColNo == 2) {
                    snpStr = string(colValue);
                    ids.push_back(snpStr);
                }
                else if (vcfColNo == 3) {
                    refStr = string(colValue);
                    refs.push_back(refStr);
                }
                else if (vcfColNo == 4) {
                    altStr = string(colValue);
                    alts.push_back(altStr);
                }
                else if (vcfColNo == 8) {
                    gtyStr = string(colValue);
                }
                else if (vcfColNo > 8)  {
                    string gtStr = string(colValue);
                    snpGts.push_back(gtStr);
                }

                valPos = 0;
                colValue[0] = 0;
                vcfColNo++;
            }
            else {
                if (buffer[buffPos] != '\r' && valPos < WORDLEN-1) {
                    colValue[valPos] = buffer[buffPos];
                    valPos++;
                }
                else {
                    colValue[valPos] = 0;
                }
            }

            if (buffer[buffPos] == '\n') {
                vcfColNo = 0;
                lineNo++;

                bool isGt = false;
                if (gtyStr.length() > 1 && gtyStr[0] == 'G' && gtyStr[1] == 'T') isGt = true;
                int numCols = snpGts.size();

                if (strcmp(chrStr.c_str(), "#CHROM") == 0) {
                    if (numSamples > 0) {
                        assert(numCols == numSamples);
                    }
                    else {
                        if (numCols > 0) {
                            for (int smpNo = 0; smpNo < numCols; smpNo++) vcfSamples.push_back(snpGts[smpNo]);
                            numSamples = numCols;
                            numGenoChars = (numSamples - 1) / 4 + 1; // Save 4 genotypes to one char
                            cout << "\tVcf file has " << numSamples << " samples\n";
                        }
                        else {
                            cout << "\nERROR: vcf file " << vcfFile << " doesn't include samples!\n";
                            fileDone = true;
                            return false;
                        }
                    }

                    hasHeadRow = true;
                    snpGts.clear();
                }
                else if (chrStr[0] && chrStr[0] != '#') {
                    if (!hasHeadRow) {
                        cout << "\nERROR: didn't find #CHROM row in vcf file\n";
                        return false;
                    }

                    if (numCols != numSamples) {
                        cout << "\nERROR at line #" << lineNo << ": chr " << chrStr
                             << ", pos " << posStr << ", snp " << snpStr
                             << ". #genotypes (" << numCols << ") is different from #samples (" << numSamples << ").\n";
                        return false;
                    }

                    int chr = GetChromosomeFromString(chrStr.c_str());
                    int rsNum = GetRsNumFromString(snpStr.c_str());
                    int pos = 0;
                    try { pos = stoi(posStr); }
                    catch (exception &err) { pos = 0; }

                    const char* gRef = refStr.c_str();
                    const char* gAlt = altStr.c_str();
                    
                    int rsSnpId = ancSnps->FindSnpIdGivenRs(rsNum);
                    int gb37SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 37, gRef, gAlt);
                    int gb38SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 38, gRef, gAlt);

                    if (isGt && (rsSnpId > -1 || gb37SnpId > -1 || gb38SnpId > -1)) {
                        putativeAncSnps++;
                        if (rsSnpId > -1)   numRsIdAncSnps++;
                        if (gb37SnpId > -1) numGb37AncSnps++;
                        if (gb38SnpId > -1) numGb38AncSnps++;

                        vcfAncSnpChrs.push_back(chr);
                        vcfAncSnpPoss.push_back(pos);
                        vcfAncSnpSnps.push_back(snpStr);
                        vcfAncSnpRefs.push_back(refStr);
                        vcfAncSnpAlts.push_back(altStr);

                        vcfRsIdAncSnpIds.push_back(rsSnpId);
                        vcfGb37AncSnpIds.push_back(gb37SnpId);
                        vcfGb38AncSnpIds.push_back(gb38SnpId);
                        
                        vector<unsigned char> gtVals;
                        for (int charNo = 0; charNo < numGenoChars; charNo++) {
                            unsigned char charGeno = 0;
                          
                            for (int charGenoNo = 0; charGenoNo < 4; charGenoNo++) {
                                int smpNo = charNo * 4 + charGenoNo;
                                if (smpNo < numSamples) {
                                    string gtStr = snpGts[smpNo];
                                    int gtLen = gtStr.length();
                                    unsigned char numAlts = 3;
                                    if (gtLen > 2 && (gtStr[1] == '|' || gtStr[1] == '/')) {
                                        char refGeno = gtStr[0];
                                        char altGeno = gtStr[2];
                                        if (refGeno == '0' && altGeno == '0') {
                                            numAlts = 0;
                                        }
                                        else if (refGeno == '0' && altGeno == '1') {
                                            numAlts = 1;
                                        }
                                        else if (refGeno == '1' && altGeno == '0') {
                                            numAlts = 1;
                                        }
                                        else if (refGeno == '1' && altGeno == '1') {
                                            numAlts = 2;
                                        }
                                    }
                                    charGeno += numAlts * charBaseInts[charGenoNo];
                                }
                            }
                            
                            gtVals.push_back(charGeno);
                        }

                        vcfAncSnpGtVals.push_back(gtVals);
                    }

                    numVcfSnps++;
                    snpGts.clear();
                }

                if (lineNo % 1000000 == 0) {
                    cout << "\tChecked " << lineNo << " lines. Found " << putativeAncSnps << " lines with ancestry SNPs\n";
                }
            }

            buffPos++;
        }

        buffNo++;
    }

    lineNo++;

    cout << "Done. Checked " << lineNo << " lines. Found " << putativeAncSnps << " lines with ancestry SNPs\n";
    gzclose (file);
    totVcfSnps += numVcfSnps;

    return true;
}

void VcfSampleAncestrySnpGeno::RecodeSnpGenotypes()
{
    // Rs ID, GB37, or GB38, use whichever returns the most ancestry SNPs to find these SNPs
    ancSnpType = AncestrySnpType::RSID;
    int maxVcfAncSnps = numRsIdAncSnps;

    if (numGb37AncSnps > maxVcfAncSnps) {
        ancSnpType = AncestrySnpType::GB37;
        maxVcfAncSnps = numGb37AncSnps;
    }

    if (numGb38AncSnps > maxVcfAncSnps) {
        ancSnpType = AncestrySnpType::GB38;
        maxVcfAncSnps = numGb38AncSnps;
    }

    int saveSnpNo = 0; // Putative SNPs saved after vcf file was read
    int ancSnpNo = 0;  // Final list of ancestry SNPs to be used for ancestry inference

    unsigned int charBaseInts[4] = {64, 16, 4, 1};
    unsigned int baseNums[8];
    unsigned int base = 1;
    for (int i = 0; i < 8; i++) {
        baseNums[i] = base;
        base = base << 1;
    }
    
    for (saveSnpNo = 0; saveSnpNo < putativeAncSnps; saveSnpNo++) {
        int ancSnpId = -1;
        if      (ancSnpType == AncestrySnpType::RSID) ancSnpId = vcfRsIdAncSnpIds[saveSnpNo];
        else if (ancSnpType == AncestrySnpType::GB37) ancSnpId = vcfGb37AncSnpIds[saveSnpNo];
        else if (ancSnpType == AncestrySnpType::GB38) ancSnpId = vcfGb38AncSnpIds[saveSnpNo];

        if (ancSnpId > -1) {
            char eRef = ancSnps->snps[ancSnpId].ref;
            char eAlt = ancSnps->snps[ancSnpId].alt;

            string vcfRef = vcfAncSnpRefs[saveSnpNo];
            string vcfAlt = vcfAncSnpAlts[saveSnpNo];

            int expRefIdx = -1;
            int expAltIdx = -1;

            CompareAncestrySnpAlleles(vcfRef, vcfAlt, eRef, eAlt, &expRefIdx, &expAltIdx);

            if (expRefIdx > -1 && expAltIdx > -1) {
                unsigned char* smpGenos = new unsigned char[numGenoChars];
                vector<unsigned char> gtVal = vcfAncSnpGtVals[saveSnpNo];

                for (int charNo = 0; charNo < numGenoChars; charNo++) {
                    int charGtVal = int(gtVal[charNo]);
                    int charGeno = 0;
                    
                    for (int cGenoNo = 0; cGenoNo < 4; cGenoNo++) {
                        int smpNo = charNo * 4 + cGenoNo;
                        if (smpNo < numSamples) {
                            int numVcfAlts = charGtVal / charBaseInts[cGenoNo];
                            charGtVal -= numVcfAlts * charBaseInts[cGenoNo];

                            char geno = 3;                            

                            if (expAltIdx == 1) {
                                geno = numVcfAlts; 
                            }
                            else if (expAltIdx == 0) {
                                geno = numVcfAlts; 
                                if (numVcfAlts < 3) {
                                    geno = 2 - numVcfAlts;
                                }
                            }

                            charGeno += geno * charBaseInts[cGenoNo];
                        }
                    }
                    
                    smpGenos[charNo] = charGeno;
                }
                
                vcfAncSnpIds.push_back(ancSnpId);
                vcfAncSnpCodedGenos.push_back(smpGenos);
                ancSnpNo++;
            }
        }
    }

    DeleteAncSnpGtValues();
}

void VcfSampleAncestrySnpGeno::CompareAncestrySnpAlleles(const string refStr, const string altsStr,
                                    const char eRef, const char eAlt, int* expRefIdx, int* expAltIdx)
{
    // Let index of the ref alleles in vcf be 0, and indices of alts be 1, 2, 3, ...
	  // Check them to find which one is the expected ref and and which is alt for the Ancestry SNP
    *expRefIdx = -1;
    *expAltIdx = -1;

    int refLen = refStr.length();
    if (refLen == 1) {
        char ref = refStr[0];
        char fRef = FlipAllele(ref);  // Flipped ref allele

        if (ref == eRef || fRef == eRef) {
            *expRefIdx = 0;
        }
        else if (ref == eAlt || fRef == eAlt) {
            *expAltIdx = 0;
        }

        // If ref seems to be flipped, then all the alts are flipped
        bool flip = false;
        if (fRef == eRef || fRef == eAlt) {
            flip = true;
        }

        if (*expRefIdx > -1 || *expAltIdx > -1) {
            vector<string> altWords = SplitString(altsStr, ",");
            int numAlts = altWords.size();

            for (int altNo = 0; altNo < numAlts; altNo++) {
                string altWord = altWords[altNo];
                int alleleIdx = altNo + 1;  // allele index starts from 0 = ref, then 1 = first alt, ...
	              if (altWord.length() == 1) {
                    char alt = altWord[0];
                    if (flip)  alt = FlipAllele(alt);

                    if      (alt == eRef && *expRefIdx < 0) *expRefIdx = alleleIdx;
                    else if (alt == eAlt && *expAltIdx < 0) *expAltIdx = alleleIdx;
                }
            }
        } // end if (*expRefIdx ...)
    }
}

int VcfSampleAncestrySnpGeno::RecodeGenotypeGivenIntegers(const int expRefIdx, const int expAltIdx, const int g1Num, const int g2Num)
{
    int genoInt = 3; // number of alt alleles, valid counts are 0, 1, 2
    
    // Genotype is valid only if both alleles are valid, i.e., same as one of the expected alleles
    int numValidGenos = 0;
    int numAlts = 0;
    
    if (g1Num == expRefIdx || g1Num == expAltIdx) {
        numValidGenos++;
        if (g1Num == expAltIdx) numAlts++;
    }
    
    if (g2Num == expRefIdx || g2Num == expAltIdx) {
        numValidGenos++;
        if (g2Num == expAltIdx) numAlts++;
    }
    
    if (numValidGenos == 2) genoInt = numAlts;
    
    return genoInt;
}

void VcfSampleAncestrySnpGeno::ShowSummary()
{
    cout << "\nTotal " << totAncSnps << " ancestry SNPs used by GrafPop\n";

    cout << "Number of samples found in the vcf file: " << numSamples << "\n";
    cout << "Total " << putativeAncSnps << " ancestry SNPs found from "	<< totVcfSnps << " SNPs\n";

    cout << "\n#RSID Ancs: " << numRsIdAncSnps << "\n"
    << "#GB37 Ancs: " << numGb37AncSnps << "\n"
    << "#GB38 Ancs: " << numGb38AncSnps << "\n";
}
