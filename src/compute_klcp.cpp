#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "ExactLCPk.hpp"
#include "NaiveLCPk.hpp"
#include "HeuristicLCPk.hpp"

void print_values(const ReadsDB& rdb,
                const ivec_t lcpKXY, const unsigned& k,
                std::ostream& lfs, const char* key){

    lfs << "\"" << key << "\"      : {" << std::endl;
    lfs << "   \"" << rdb.getReadNameById(0)<< "\": [" << std::endl;
    for(unsigned i = 0; i < lcpKXY.size(); i++)
        lfs << "        "
            << "[" << i
            << ", " << lcpKXY[i]
            << ((i == lcpKXY.size() - 1) ? "]" : "],\n");
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
    cfg.ofs << "{" << std::endl;
    ivec_t llrk = lxy.getkLCP()[0][1];
    print_values(rdb, llrk, cfg.kv, cfg.ofs, "llrk");

    unsigned size = s.size();
    ivec_t workspace(size, 0);
    // calculate LSUS
    for(unsigned i = 0; i < size; ++i) {
        if(llrk[i] + i == size) {
            llrk[i] = -1;
        } else {
            llrk[i] = llrk[i] + i;
        }
    }
    ivec_t lsus = llrk;
    print_values(rdb, lsus, cfg.kv, cfg.ofs, "lsusk");

    // calculate SLS

    
    // calculate SUS
    
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
