//
// Created by shifeng on 6/10/21.
//
#include "Memusage.h"

void process_mem_usage() {
double vm_usage     = 0.0;
double resident_set = 0.0;

// the two fields we want
unsigned long vsize;
long rss;
{
std::string ignore;
std::ifstream ifs("/proc/self/stat", std::ios_base::in);
ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
>> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
>> ignore >> ignore >> vsize >> rss;
}

long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
vm_usage = vsize / 1024.0;
resident_set = rss * page_size_kb;

std::cout << "       VM (MB): " << int(vm_usage/1024) << "       RSS (MB): " << int(resident_set/1024) << std::endl;
}
