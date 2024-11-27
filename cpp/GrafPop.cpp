#include "GrafPop.h"

SampleGenoAncestry *smpGenoAnc = NULL;

int main(int argc, char* argv[])
{
    string usage = "Usage: grafanc <Binary PLINK set or VCF file> <output file>\n";

    string disclaimer =
    "\n *==========================================================================="
    "\n *  GrafAnc: Software to infer subject ancestry from genotypes quickly"
    "\n *  Yumi (Jimmy) Jin, PhD"
    "\n *  Jimmy.Jin@Pennmedicine.upenn.edu"
    "\n *  10/08/2024"
    "\n *"
    "\n *===========================================================================";

    if (argc < 3) {
        cout << usage << "\n";
        exit(0);
    }
    
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    string genoDs, outputFile;
    genoDs = argv[1];
    outputFile = argv[2];
    
    parameters params = GetParameters(argc, argv);
    
    string fileBase = "";
    cout << "Checking file " << genoDs << "\n";
    
    GenoDatasetType fileType = CheckGenoDataFile(genoDs, &fileBase);

    if (fileType == GenoDatasetType::NOT_EXISTS) {
        cout << "\nERROR: Genotype file " << genoDs << " doesn't exist!\n\n";
    	  return 0;
    }
    else if (fileType == GenoDatasetType::IS_PLINK_GZ) {
        cout << "\nERROR: PLINK set " << genoDs << " is zipped. Please unzip it.\n\n";
        return 0;
    }
    else if (fileType == GenoDatasetType::IS_OTHER) {
        cout << "\nERROR: Genotype file " << genoDs << " should be a binary PLINK set or vcf or vcf.gz file..\n\n";
        return 0;
    }

    string ancSnpFile = FindFile("AncSnpPopAFs.txt");
    if (ancSnpFile == "") {
        cout << "ERROR: didn't find file AncSnpPopAFs.txt. Please put the file under 'data' directory.\n\n";
        return 0;
    }

    string refSubPopFile = ancSnpFile;
    string nomSubPopFile = "";
    
    AncestrySnps *ancSnps = new AncestrySnps();
    int numSnpsInAncFile = ancSnps->ReadAncestrySnpsFromFile(ancSnpFile);

    int numSnpsInRefPopFile = ancSnps->ReadRefSubPopSnpsFromFile(refSubPopFile);
    cout << "Read " << numSnpsInRefPopFile << " SNPs from " << refSubPopFile << "\n"; 
    
    int totAncSnps = ancSnps->GetNumAncestrySnps();
    int minAncSnps = 100;

    int numThreads = thread::hardware_concurrency();
    if (params.numThreads && params.numThreads < numThreads) numThreads = params.numThreads;
    
    smpGenoAnc = new SampleGenoAncestry(ancSnps, minAncSnps);
    smpGenoAnc->CalculateSubPopGd0Values();

    //// Get number of samples from PLINK fam file or vcf file
    string bedFile = "";
    string bimFile = "";
    string famFile = "";
    
    FamFileSamples *famSmps = NULL;
    BedFileSnpGeno *bedGeno = NULL;
    BimFileAncestrySnps *bimSnps = NULL;
    int numBimAncSnps = 0;
    
    int numDsSamples = 0;
    string fileTypeStr = "";
    
    if (fileType == GenoDatasetType::IS_VCF || fileType == GenoDatasetType::IS_VCF_GZ) {
        VcfGrafAncSnpGeno* vcfHead = new VcfGrafAncSnpGeno(genoDs, ancSnps);
        bool headRead = vcfHead->ReadHeaderFromFile();
        numDsSamples = vcfHead->GetTotalVcfSamples();
        fileTypeStr = "VCF";
        
        delete vcfHead;
    }
    else if (fileType == GenoDatasetType::IS_PLINK) {
        bedFile = fileBase + ".bed";
        bimFile = fileBase + ".bim";
        famFile = fileBase + ".fam";
      
        if ( !FileExists(bedFile.c_str()) ||
             !FileExists(bimFile.c_str()) ||
             !FileExists(famFile.c_str())    ) {
             if (!FileExists(bedFile.c_str())) cout << "\nERROR: didn't find " << bedFile << "\n";
             if (!FileExists(bimFile.c_str())) cout << "\nERROR: didn't find " << bimFile << "\n";
             if (!FileExists(famFile.c_str())) cout << "\nERROR: didn't find " << famFile << "\n";
             cout << "\n";
             
             return 0;
        }
        
        famSmps = new FamFileSamples(famFile);
        famSmps->ShowSummary();
        
        bimSnps = new BimFileAncestrySnps(totAncSnps);
        bimSnps->ReadAncestrySnpsFromFile(bimFile, ancSnps);
        
        numBimAncSnps = bimSnps->GetNumBimAncestrySnps();
        
        bimSnps->ShowSummary();
        
        numDsSamples = famSmps->GetNumFamSamples();
        fileTypeStr = "PLINK";
    }
    
    cout << "Number of samples in " << fileTypeStr << " dataset: " << numDsSamples << "\n";
    if (fileTypeStr == "PLINK") cout << "File is PLINK\n";    

    //// Check available memory
    int availMem = GetAvailableMemoryInMb();
    
    int maxMem = 8000;
    if (params.maxMemMb) maxMem = params.maxMemMb;
    
    int totAllocMem = GetAllocatableMemoryInMb();
    cout << "Available memory: " << totAllocMem << " MiB\n";
    if (totAllocMem > maxMem) {
        cout << "Maximum "  << maxMem << " MiB will be used.\n";
        totAllocMem = maxMem;
    }
    
    float genosPerByte = fileTypeStr == "PLINK" ? 2 : 1;

    //// Determine how many rounds are needed for analyzing all samples
    int memNeeded = int( (numDsSamples/1000) * (totAncSnps/1000) / genosPerByte );

    int numSmpBlocks = 1;
    int smpBlockSize = numDsSamples;

    cout << "Need memory: " << memNeeded << " MiB\n";
    if (memNeeded > totAllocMem) {    
        smpBlockSize = int(totAllocMem * genosPerByte * 1000000 / totAncSnps);
        smpBlockSize = int(smpBlockSize / 100.0) * 100; // Bump down to 100's
        if (params.blockSize && params.blockSize < smpBlockSize) smpBlockSize = params.blockSize;
        numSmpBlocks = (numDsSamples-1) / smpBlockSize + 1;
        
        cout << numSmpBlocks << " rounds are needed to analyze the " << numDsSamples << " samples. Each round analyzes "  << smpBlockSize << " samples.\n";
    }

    VcfGrafAncSnpGeno *vcfGeno = NULL;
    
    for (int round = 0; round < numSmpBlocks; round++) {
        if (numSmpBlocks > 1) cout << "\nRound No. " << round+1 << "\n";
      
        struct timeval rt1, rtm, rt2;
        gettimeofday(&rt1, NULL);
        
        bool hasVcfGeno = false;
        bool hasBedGeno = false;
        
        if (fileType == GenoDatasetType::IS_VCF || fileType == GenoDatasetType::IS_VCF_GZ) {
            vcfGeno = new VcfGrafAncSnpGeno(genoDs, ancSnps);
            hasVcfGeno = true; 
    
            bool verbose = round == 0 ? true : false;
            bool dataRead = vcfGeno->ReadDataFromFile(round*smpBlockSize, smpBlockSize, verbose);

            if (!dataRead) {
                cout << "\nFailed to read genotype data from " << genoDs << "\n\n";
                return 0;
            }
            vcfGeno->RecodeSnpGenotypes();
            
            int totVcfSnps = vcfGeno->GetNumVcfSnps();
            int numAncSnps = vcfGeno->vcfAncSnpIds.size();
            
            int numVcfSmps = vcfGeno->GetTotalVcfSamples();
            int numChkSmps = vcfGeno->GetNumSamples();
            
            if (round == 0) cout << "Total " << totVcfSnps << " SNPs in VCF file. " << numAncSnps << " SNPs are ancestry SNPs.\n"; 
            if (round == 0) cout << "Total " << numVcfSmps << " samples in file.\n";
            if (numSmpBlocks > 1) cout << "\tGenotypes read for " << numChkSmps << " samples.\n";

            if (smpGenoAnc->HasEnoughAncestrySnps(numAncSnps)) {
                smpGenoAnc->SetGenoSamples(vcfGeno->vcfSamples);
                smpGenoAnc->SetSnpGenoData(&vcfGeno->vcfAncSnpIds, &vcfGeno->vcfAncSnpCodedGenos);
            }
            else {
                cout << "\nWARNING: Ancestry inference not done due to lack of genotyped ancestry SNPs "
                 << "(at least " << minAncSnps << " ancestry SNPs are needed).\n\n";
                return 0;
            }
        }
        else if (fileType == GenoDatasetType::IS_PLINK) {
            vector<FamSample> chkSmps;
            int stSmpNo = round * smpBlockSize;
            int numChkSmps = smpBlockSize;
            if (stSmpNo + numChkSmps > numDsSamples) numChkSmps = numDsSamples - stSmpNo;
    
            for (int i = stSmpNo; i < stSmpNo+numChkSmps; i++) chkSmps.push_back(famSmps->samples[i]);
            
            smpGenoAnc->SetGenoSamples(chkSmps);
            
            int numSmps = smpGenoAnc->GetNumSamples();

            if (smpGenoAnc->HasEnoughAncestrySnps(numBimAncSnps)) {
                bedGeno = new BedFileSnpGeno(bedFile, ancSnps, bimSnps, famSmps);
                hasBedGeno = true;
    
                bedGeno->SelectFamSampleIds(stSmpNo, numChkSmps);
                bool hasErr = bedGeno->ReadGenotypesFromBedFile();
                
                if (hasErr) return 0;
                bedGeno->ShowSummary();
                
                smpGenoAnc->SetSnpGenoData(&bedGeno->ancSnpSnpIds, &bedGeno->ancSnpSmpGenos);
            }
            else {
                cout << "Ancestry inference not done due to lack of genotyped ancestry SNPs.\n\n";
                return 0;
            }
        }
        
        cout << "\tDone reading genotypes.\n";
        gettimeofday(&rtm, NULL);
        ShowTimeDiff(rt1, rtm);
        
        cout << "Launching " << numThreads << " threads to calculate ancestry scores.\n";
        smpGenoAnc->SetNumThreads(numThreads);
    
        mutex iomutex;
        vector<thread> threads(numThreads);
    
        for (unsigned i = 0; i < numThreads; ++i) {
            threads[i] = thread([&iomutex, i] {
                {
                    lock_guard<mutex> iolock(iomutex);
                }
    
      	        smpGenoAnc->SetAncestryPvalues(i);
            });
        }
    
        for (auto& t : threads) {
            t.join();
        }
    
        //// Append results to the output file
        bool isAppend = round == 0 ? false : true;
        smpGenoAnc->SaveAncestryResults(outputFile, isAppend);
    
        if (hasVcfGeno) delete vcfGeno;
        if (hasBedGeno) delete bedGeno;
        
        smpGenoAnc->ResetSamples();

        cout << "\tResults saved to " << outputFile << "\n";
        gettimeofday(&rt2, NULL);
        ShowTimeDiff(rtm, rt2);
    }
    cout << "\n";
    
    delete famSmps;
    delete bimSnps;
    
    gettimeofday(&t2, NULL);
    ShowTimeDiff(t1, t2);

    return 1;
}

string GetExecutablePath()
{
    char rawPathName[PATH_MAX];
    realpath(PROC_SELF_EXE, rawPathName);

    string exePath = string(rawPathName);
    size_t slashPos = exePath.find_last_of("/\\");
    string exeDir = exePath.substr(0, slashPos);

    return exeDir;
}

string FindFile(string filename)
{
    string fullFile = filename;

    if (FileExists(fullFile.c_str())) return fullFile;

    string exeDir = GetExecutablePath();
    fullFile = exeDir + "/data/" + filename;
    if (FileExists(fullFile.c_str())) return fullFile;

    fullFile = exeDir + "/" + filename;
    if (FileExists(fullFile.c_str())) return fullFile;

    if(const char* grafPath = getenv("GRAFPATH")) {
        string grafDir = string(grafPath);
        fullFile = grafDir + "/data/" + filename;
        if (FileExists(fullFile.c_str())) return fullFile;

        fullFile = grafDir + "/" + filename;
        if (FileExists(fullFile.c_str())) return fullFile;
    }

    return "";
}

int GetAvailableMemoryInMb()
{
    struct sysinfo memInfo;
    
    sysinfo (&memInfo);
    long long totVirtMem = memInfo.totalram;
    totVirtMem += memInfo.totalswap;
    totVirtMem *= memInfo.mem_unit;
    
    int totPhysMem = memInfo.totalram / 1000000;
    totPhysMem *= memInfo.mem_unit;
    
    return totPhysMem;
}

int GetAllocatableMemoryInMb()
{
    int mbSize = 1000000; 
    size_t blockSize = mbSize * 100;
    size_t totalAllocated = blockSize * 1000;
    
    int iter = 0;
    bool done = false;
    
    while (!done) {
        try {
            char* ptr = new char[totalAllocated];
            delete[] ptr;
            done = true;
        } 
        catch (std::bad_alloc& ba) {
            totalAllocated -= blockSize;
        }
        
        iter++;
        if (iter > 2000) done = true;
    }
  
    int mbAllocated = totalAllocated / mbSize;  
    
    return mbAllocated;
}

parameters GetParameters(int argc, char** argv)
{
    parameters param = {0, 0, 0};

    int optRes;
    int optIdx;

    static struct option options[] =
        {
            {"maxmem",   required_argument, 0, 0},
            {"block",    required_argument, 0, 0},
            {"threads",  required_argument, 0, 0},
            {0, 0, 0, 0}
        }; 
    
    
    while ((optRes = getopt_long(argc, argv, "", options, &optIdx)) != -1)
    {
        switch (optRes)
        {
            // 0 means long option.
            case 0:
            {
                if (optarg) {
                    int argVal = atoi(optarg);
                    if (strcmp(options[optIdx].name, "maxmem") == 0)  param.maxMemMb   = argVal;
                    if (strcmp(options[optIdx].name, "block") == 0)   param.blockSize  = argVal;
                    if (strcmp(options[optIdx].name, "threads") == 0) param.numThreads = argVal;
                }    
                break;
            }
        }
    }
    
    // Set minimum requied memory and block size
    if (param.maxMemMb > 0 && param.maxMemMb < 100) param.maxMemMb = 100;
    if (param.blockSize > 0 && param.blockSize < 100) param.blockSize = 100;
    param.blockSize = int(param.blockSize / 100.0) * 100; // Bump down to 100's
    
    return param;
}


