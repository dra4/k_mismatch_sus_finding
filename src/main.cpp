#include "ReadsDB.hpp"
#include "AppConfig.hpp"
#include "compute_klcp.hpp"
#include "compute_kacs.hpp"

#include <iostream>
#include <ctime>

int main(int argc, char** argv){
    AppConfig acfg(argc, argv);
    if(!acfg.validate(std::cout))
        return 0;
    clock_t startTime = clock();

    // Build the reads database and get it as a string
    ReadsDB rdb(acfg.dir, acfg.ifiles, 0);

    if(rdb.getReadsCount() == 0) {
        std::cout << "\"error\" : \"NO READS AVAILABLE\"" << std::endl;
        return 0;
    }
    // write basic config/input meta data to log file
    acfg.lfs << "{" << std::endl;
    acfg.write(acfg.lfs);
    rdb.writeMeta(acfg.lfs);

    clock_t klcpStart = clock();

    // estimate k-lcp
    if(acfg.only_lcp){
        compute_klcp(rdb, acfg);
    } else {
        compute_kacs(rdb, acfg);
    }

    // write time taken
    acfg.lfs << "\"klcp_time_secs\" : "
             << double( clock() - klcpStart ) / ((double)CLOCKS_PER_SEC )
             << "," << std::endl;
    acfg.lfs << "\"total_time_secs\" : "
             << double( clock() - startTime ) / ((double)CLOCKS_PER_SEC )
             << std::endl;
    acfg.lfs << "}" << std::endl;
    return 0;
}
