#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "ExactLCPk.hpp"
#include "NaiveLCPk.hpp"
#include "HeuristicLCPk.hpp"

void print_lcpk(const ReadsDB& rdb,
                const ivec_t lcpKXY[2][2], const unsigned& k,
                std::ostream& lfs, const char* key){
    /* const std::string& sx = rdb.getReadById(i); */
    /* const std::string& sy = rdb.getReadById(j); */

    lfs << "\"" << key << "\"      : {" << std::endl;
/* #ifdef DEBUG */
    /* lfs << "   \"meta\" : [" << i << ",\t" << j << ",\t" << sx.size() << ",\t" */
    /*     << sy.size() << ",\t" << k << "]," << std::endl; */
/* #endif */
    lfs << "   \"" << rdb.getReadNameById(0)<< "\": [" << std::endl;
    for(unsigned i = 0; i < lcpKXY[0][0].size(); i++)
        lfs << "        "
            << "[" << i
            << ", " << lcpKXY[0][0][i]
            << ", " << lcpKXY[0][1][i]
            << ((i == lcpKXY[0][0].size() - 1) ? "]" : "],\n");
    /* lfs << "]," << std::endl; */
    /* lfs << "   \"" << rdb.getReadNameById(j)<< "\": [" << std::endl; */
    /* for(unsigned i = 0; i < lcpKXY[1][0].size(); i++) */
    /*     lfs << "        " */
    /*         << "[" << i */
    /*         << ", " << lcpKXY[1][0][i] */
    /*         << ", " << lcpKXY[1][1][i] */
    /*         << ((i == lcpKXY[1][0].size() - 1) ? "]" : "],\n"); */
    lfs << "]" << std::endl;
    lfs << "  }," << std::endl;

}

template<typename LCPk>
void klcp_pair_factory(ReadsDB& rdb, AppConfig& cfg){
#ifdef DEBUG
    cfg.lfs << "\"klcp_debug\"      : [" << std::endl;
#endif
    const std::string& s = rdb.getReadById(0);

    LCPk lxy(s, s, cfg); // construct suffix array
#ifdef DEBUG
    lxy.print(cfg.lfs);
#endif
    lxy.compute();
#ifdef DEBUG
    cfg.lfs << " []]," << std::endl;
#endif
#ifndef DEBUG_KLCP
    print_lcpk(rdb, lxy.getkLCP(), 1, cfg.lfs, "llrk");
#endif
    cfg.ofs << "{" << std::endl;
    print_lcpk(rdb, lxy.getkLCP(), cfg.kv, cfg.ofs, "llrk");
    cfg.ofs << "  \"end\": []" << std::endl
            << "}" << std::endl;

}

void compute_klcp(ReadsDB& rdb, AppConfig& cfg){
    // TODO
    unsigned nReads = rdb.getReadsCount();
    assert(nReads >= 1);

    klcp_pair_factory<ExactLCPk>(rdb, cfg);
    /* for(unsigned i  = 0; i < nReads; i++){ */
    /*     for(unsigned j = i + 1; j < nReads; j++){ */
    /*         if(cfg.method == 1) */
    /*             klcp_pair_factory<NaiveLCPk>(i, j, rdb, cfg); */
    /*         else if(cfg.method == 2) */
    /*             klcp_pair_factory<HeuristicLCPk>(i, j, rdb, cfg); */
    /*         else */
    /*             klcp_pair_factory<ExactLCPk>(i, j, rdb, cfg); */
    /*         break; */
    /*     } */
    /*     break; */
    /* } */
}
