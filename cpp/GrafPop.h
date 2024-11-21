using namespace std;

#ifndef GRAFPOP_H
#define GRAFPOP_H

#include "Util.h"
#include "AncestrySnps.h"
#include "VcfGrafAncSnpGeno.h"
#include "FamFileSamples.h"
#include "BimFileAncestrySnps.h"
#include "BedFileSnpGeno.h"
#include "SampleGenoDist.h"
#include "SampleGenoAncestry.h"
#include <thread>
#include <mutex>

#include "htslib/hts.h"
#include "htslib/vcf.h"

#include "sys/types.h"
#include "sys/sysinfo.h"

#if defined(__sun)
#define PROC_SELF_EXE "/proc/self/path/a.out"
#else
#define PROC_SELF_EXE "/proc/self/exe"
#endif

string GetExecutablePath(void);
string FindFile(string);

int GetAvailableMemoryInMb();
int GetAllocatableMemoryInMb();

#endif
