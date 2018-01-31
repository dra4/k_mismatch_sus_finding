#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctime>

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

void print_a_b(const ivec_t a, const ivec_t b) {
    std::ofstream ab_out;
    ab_out.open("alfred_a_b.txt", std::ofstream::out);
    for(int i = 0; i < a.size(); ++i) {
        ab_out << a[i] << ", " << b[i] << "\n";
    }
    ab_out.close();
}

void print_b(const ivec_t b) {
    std::ofstream b_out;
    b_out.open("alfred_b.txt", std::ofstream::out);
    for(int i = 0; i < b.size(); ++i) {
        b_out << b[i] << "\n";
    }
    b_out.close();
}

template<typename LCPk>
void klcp_pair_factory(ReadsDB& rdb, AppConfig& cfg){
#ifdef DEBUG
    cfg.lfs << "\"klcp_debug\"      : [" << std::endl;
#endif
    const std::string& s = rdb.getReadById(0);
    clock_t startTime = clock();

    LCPk lxy(s, s, cfg); // construct suffix array
#ifdef DEBUG
    lxy.print(cfg.lfs);
#endif
    lxy.compute();
#ifdef DEBUG
    cfg.lfs << " []]," << std::endl;
#endif
    cfg.ofs << "{" << std::endl;
    /* std::cout << "Getting LLRk...\n"; */
    ivec_t &llrk = lxy.getkLCP()[0][1];
    /* print_values(rdb, llrk, cfg.kv, cfg.ofs, "llrk"); */

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
    ivec_t &lsusk = llrk;
    /* print_b(lsusk); */
    /* print_values(rdb, lsusk, cfg.kv, cfg.ofs, "lsusk"); */

    // calculate SLS
    /* std::cout << "calculating SLSk...\n"; */
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
    /* print_values(rdb, workspace, cfg.kv, cfg.ofs, "pred"); */

    for(int32_t i = 0; i <= r; ++i) {
        if(workspace[i] == NIL) {
            workspace[i] = i;
        } else {
            int32_t val = lsusk[workspace[i]] + 1;
            workspace[i] = val > i ? val : i;
        }
    }
    /* print_values(rdb, workspace, cfg.kv, cfg.ofs, "t"); */

    for(int32_t i = r; i >= 0; --i) {
        if(workspace[workspace[i]] == NIL) {
            workspace[workspace[i]] = i;
        }
        if(i < workspace[i]) {
            workspace[i] = NIL;
        }
    }

    /* print_values(rdb, workspace, cfg.kv, cfg.ofs, "t_1"); */

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
    /* print_b(workspace); */
    /* print_values(rdb, workspace, cfg.kv, cfg.ofs, "slsk"); */

    
    // calculate SUS
    /* std::cout << "calculating SUSk...\n"; */
    int32_t z = NIL;
    for(int32_t i = size - 1; i >= 0; --i) {
        if(lsusk[i] != NIL) {
            z = lsusk[i];
            break;
        }
    }

    if(z < size - 1) {
        for(int32_t i = z + 1; i < size; ++i) {
            lsusk[i] = NIL;
        }
    }

    for(int32_t i = z; i >= 0; --i) {
        lsusk[i] = lsusk[workspace[i]];
    }

    for(int32_t i = 1; i < size; ++i) {
        if(workspace[i] == NIL && lsusk[i] == NIL) {
            workspace[i] = workspace[i - 1];
            lsusk[i] = lsusk[i-1] + 1;
        } else if(lsusk[i-1] == i - 1 && ((lsusk[i-1] - workspace[i-1] + 2) < (lsusk[i] - workspace[i] + 1))) {
            workspace[i] = workspace[i-1];
            lsusk[i] = lsusk[i-1] + 1;
        }
    }

    std::cout << "Seconds: " << double( clock() - startTime ) / ((double)CLOCKS_PER_SEC ) << "\n";
    print_values(rdb, workspace, cfg.kv, cfg.ofs, "suska");
    print_values(rdb, lsusk, cfg.kv, cfg.ofs, "suskb");
    print_a_b(workspace, lsusk);



    
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
