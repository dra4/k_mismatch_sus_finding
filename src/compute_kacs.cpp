#include <cmath>
#include <iomanip>
#include "util.hpp"
#include "NaiveLCPk.hpp"
#include "ExactLCPk.hpp"
#include "HeuristicLCPk.hpp"

double correct_term(double xLen){
    if(xLen > 0)
        return 2 * (std::log10(xLen)) / xLen;
    return 0;
}

void compute_kacs(const ivec_t lcpkXY[2][2],
                  double& acsxy, double& acsyx){
    acsxy = 0;
    for(unsigned i = 0; i < lcpkXY[0][1].size(); i++){
        acsxy += lcpkXY[0][1][i];
    }
    acsxy = acsxy / (double)lcpkXY[0][1].size();

    acsyx = 0;
    for(unsigned i = 0; i < lcpkXY[1][1].size(); i++){
        acsyx += lcpkXY[1][1][i];
    }
    acsyx = acsyx / (double)lcpkXY[1][1].size();
}

void compute_kdist(const double& xlen, const double& ylen,
                   const double& acsxy, const double& acsyx,
                   double& dxy){
    double dpxy, dpyx;
    dpxy = std::log10(ylen) / acsxy;
    dpxy -= correct_term(xlen);
    dpyx = std::log10(xlen) / acsyx;
    dpyx -= correct_term(ylen);
    dxy = (dpxy + dpyx)/2.0;
}


void compute_kdist(const ivec_t lcpkXY[2][2], double& dxy){
    double acsxy, acsyx;
    compute_kacs(lcpkXY, acsxy, acsyx);
    double xlen = (double)lcpkXY[0][1].size();
    double ylen = (double)lcpkXY[1][1].size();
    compute_kdist(xlen, ylen, acsxy, acsyx, dxy);
}

void write_dmat(std::vector<std::vector<double>>& dmat,
                ReadsDB& rdb, AppConfig& cfg)
{
    unsigned nReads = dmat.size();
    cfg.ofs << nReads << std::endl;
    for(unsigned i=0; i < nReads; i++){
        // removeExtension(basename(rdb.getReadNameById(i)));
        std::string orgName = rdb.getReadNameById(i);
        if(orgName.size() > 10)
            orgName = orgName.substr(0, 10);
        cfg.ofs << std::setw(14) << std::left << orgName;
        for(unsigned j = 0; j < nReads; j++){
            double kdxy = dmat[i][j];
            cfg.ofs << kdxy << ((j == nReads - 1) ? "" : " ");
        }
        cfg.ofs << std::endl;
    }

}

template<typename LCPk>
void kacs_factory(ReadsDB& rdb, AppConfig& cfg){
    unsigned nReads = rdb.getReadsCount();
    assert(nReads >= 2);
    assert(nReads == rdb.getFilesCount());
    std::vector<std::vector<double>> dmat;
    dmat.resize(nReads);
    for(unsigned i = 0;i < nReads; i++)
        dmat[i].resize(nReads);

    cfg.lfs << "\"rdata\": [" << std::endl;

    for(int i  = nReads - 1; i >= 0; i--){
        for(int j = i + 1; j < (int)nReads; j++){
            const std::string& sx = rdb.getReadById(i);
            const std::string& sy = rdb.getReadById(j);
            double xlen = sx.size(), ylen = sy.size();
            double acsxy, acsyx, kdxy;
            LCPk lxy(sx, sy, cfg); // construct suffix array
            lxy.compute();
            compute_kacs(lxy.getkLCP(), acsxy, acsyx);
            compute_kdist(xlen, ylen, acsxy, acsyx, kdxy);
            // cfg.ofs << i << "\t" << j << "\t"
            //         << sx.size() << "\t" << sy.size() << "\t"
            //         << acsxy << "\t" << acsyx << "\t"
            //         << kdxy << std::endl;
            dmat[i][j] = dmat[j][i] = kdxy;
        }
    }

    cfg.lfs << "  []]," << std::endl;
    write_dmat(dmat, rdb, cfg);
}

void compute_kacs(ReadsDB& rdb, AppConfig& cfg)
{
    // estimate k-acs
    if(cfg.method == 1)
        kacs_factory<NaiveLCPk>(rdb, cfg);
    else if(cfg.method == 2)
        kacs_factory<HeuristicLCPk>(rdb, cfg);
    else
        kacs_factory<ExactLCPk>(rdb, cfg);
}
