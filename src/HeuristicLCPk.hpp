#ifndef HEURISTIC_LCPK_H
#define HEURISTIC_LCPK_H
#include "ExactLCPk.hpp"

class HeuristicLCPk{
private:
    ExactLCPk m_eLCPk;
    AppConfig& m_aCfg;
    const std::string& m_sx;
    const std::string& m_sy;
    ivec_t m_klcpXY[2][2];
    int32_t m_strLengths[2];
    int32_t m_shiftPos[2][2];
    int m_kextend;
    int32_t rmq(int tidx, const int32_t& t1, const int32_t& t2);
    void extend(const std::string& sa, const std::string& sb, int tidx);
    void extendRMQ(const std::string& sa, const std::string& sb, int tidx);
    void extendCrawl(const std::string& sa, const std::string& sb, int tidx);
    void extendRMQ();
    void extendCrawl();
    void computeBasisTest(int kv);
    void computeBasis();
public:
    HeuristicLCPk(const std::string& x, const std::string& y,
                  AppConfig& cfg);
    void computeCrawl();
    void compute();
    auto getkLCP() -> const ivec_t (&)[2][2] {
        return m_klcpXY;
    }
    void computeCrawlTest(int kv, int ext);
    void computeTest(int kv, int ext);
    void print(std::ostream& ){}
};

#endif /* HUERISTICLCPK_H */
