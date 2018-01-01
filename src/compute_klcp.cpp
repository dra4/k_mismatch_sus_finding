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
    const int32_t NIL = -1;
    // calculate LSUS
    for(unsigned i = 0; i < size; ++i) {
        if(llrk[i] + i == size) {
            llrk[i] = NIL;
        } else {
            llrk[i] = llrk[i] + i;
        }
    }
    ivec_t lsusk = llrk;
    print_values(rdb, lsusk, cfg.kv, cfg.ofs, "lsusk");

    // calculate SLS
    int32_t r = size - 1;
    for(;r >= 0; --r) {
        if(lsusk[r] != NIL) {
            break;
        }
    }

    if(r < size - 1) {
        for(int32_t i = r + 1; i < size; ++i) {
            workspace[i] = NIL;
        }
    }

    workspace[0] = NIL;

    for(int32_t i = 1; i <= r; ++i) {
        int32_t li = lsusk[i] - i;
        int32_t j = i - 1;
        while(workspace[j] != NIL && lsusk[j] - j >= li) {
            j = workspace[j];
        }
        if(lsusk[j] - j < li) {
            workspace[i] = j;
        } else {
            workspace[i] = NIL;
        }
    }
    print_values(rdb, workspace, cfg.kv, cfg.ofs, "pred");

    for(int32_t i = 0; i <= r; ++i) {
        if(workspace[i] == NIL) {
            workspace[i] = i;
        } else {
            int32_t val = lsusk[workspace[i]] + 1;
            workspace[i] = val > i ? val : i;
        }
    }
    print_values(rdb, workspace, cfg.kv, cfg.ofs, "t");

    for(int32_t i = r; i >= 0; --i) {
        if(workspace[workspace[i]] == NIL) {
            workspace[workspace[i]] = i;
        }
        if(i < workspace[i]) {
            workspace[i] = NIL;
        }
    }

    print_values(rdb, workspace, cfg.kv, cfg.ofs, "t_1");

    workspace[0] = 0;
    for(int32_t i = 1; i <= lsusk[r]; ++i) {
        int32_t prevSls = workspace[i-1];
        int32_t curT_1 = workspace[i];
        if(prevSls > curT_1) {
            workspace[i] = prevSls;
        } else {
            workspace[i] = curT_1;
        }
    }
    print_values(rdb, workspace, cfg.kv, cfg.ofs, "slsk");

    
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
