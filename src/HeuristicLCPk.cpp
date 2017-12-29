#include "HeuristicLCPk.hpp"

HeuristicLCPk::HeuristicLCPk(const std::string& x, const std::string& y,
                             AppConfig& cfg): m_eLCPk(x, y, cfg), m_aCfg(cfg),
                                              m_sx(x), m_sy(y) {
    m_kextend = m_aCfg.extend - m_aCfg.kv;
    m_strLengths[0] = m_sx.size(); m_strLengths[1] = m_sy.size();
    m_shiftPos[0][0] = 0; m_shiftPos[0][1] = m_strLengths[0] + 1;
    m_shiftPos[1][0] = m_strLengths[0] + 1; m_shiftPos[1][1] = 0;
}


void HeuristicLCPk::extendCrawl(const std::string& sa, const std::string& sb,
                                int tidx){
    if(m_kextend <= 0)
        return;
    //std::cout << sa << std::endl;
    //std::cout << sb << std::endl;
    for(unsigned astart = 0; astart < sa.size(); astart++){
        unsigned bstart = (unsigned)m_klcpXY[tidx][0][astart];
        unsigned howfar = (unsigned)m_klcpXY[tidx][1][astart];
        int kdx = 0;
        while((bstart + howfar < sb.size()) &&
              (astart + howfar < sa.size())){
            if(sa[astart + howfar] != sb[bstart + howfar])
                kdx += 1;
            if(m_kextend < kdx)
                break;
            howfar++;
        }
        if(m_klcpXY[tidx][1][astart] < (int32_t)howfar){
            m_klcpXY[tidx][1][astart] = (int32_t)howfar;
        }
    }
}

int32_t HeuristicLCPk::rmq(int tidx, const int32_t& t1, const int32_t& t2){
    return m_eLCPk.sufRangeMinLCP(t1 + m_shiftPos[tidx][0],
                                  t2 + m_shiftPos[tidx][1]);
}


void HeuristicLCPk::extendRMQ(const std::string& sa, const std::string& sb,
                              int tidx){
    if(m_kextend <= 0)
        return;
    //std::cout << sa << std::endl;
    //std::cout << sb << std::endl;
    for(int32_t astart = 0; astart < (int32_t)sa.size(); astart++){
        int32_t bstart = m_klcpXY[tidx][0][astart];
        int32_t howfar = m_klcpXY[tidx][1][astart] - 1; // start back a bit
#ifdef DEBUG
        m_aCfg.lfs << " [" << astart << ",\t" << bstart << ",\t"
                   << howfar << ",\t";
#endif
        for(int kdx = 0; kdx <= m_kextend ; kdx++){
            if((bstart + howfar >= (int32_t)sb.size()) ||
               (astart + howfar >= (int32_t)sa.size()))
                break;
            howfar += 1;
            int32_t hx = rmq(tidx, astart + howfar, bstart + howfar);
            howfar += hx;
#ifdef DEBUG
            m_aCfg.lfs << howfar << ",\t" << hx << ",\t";
#endif
        }

        if(m_klcpXY[tidx][1][astart] < (int32_t)howfar){
            m_klcpXY[tidx][1][astart] = (int32_t)howfar;
        }
#ifdef DEBUG
        m_aCfg.lfs << howfar << "]," << std::endl;
#endif
    }
#ifdef DEBUG
    m_aCfg.lfs << std::endl << std::endl;
#endif

}

void HeuristicLCPk::computeBasis(){
    // compute
    m_eLCPk.compute();
    // sawp
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
            std::swap(m_klcpXY[i][j], m_eLCPk.m_klcpXY[i][j]);
}

void HeuristicLCPk::computeBasisTest(int kv){
    // compute
    m_eLCPk.computeTest(kv);
    // sawp
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
            std::swap(m_klcpXY[i][j], m_eLCPk.m_klcpXY[i][j]);
}

void HeuristicLCPk::extendRMQ(){
    extendRMQ(m_sx, m_sy, 0);
    extendRMQ(m_sy, m_sx, 1);
}

void HeuristicLCPk::extendCrawl(){
    extendCrawl(m_sx, m_sy, 0);
    extendCrawl(m_sy, m_sx, 1);
}

void HeuristicLCPk::compute(){
    computeBasis();
    extendRMQ();
}

void HeuristicLCPk::computeCrawl(){
    computeBasis();
    extendCrawl();
}

void HeuristicLCPk::computeTest(int kv, int ext){
    m_kextend = ext - kv;
    computeBasisTest(kv);
    extendRMQ(); // extend
}

void HeuristicLCPk::computeCrawlTest(int kv, int ext){
    m_kextend = ext - kv;
    computeBasisTest(kv);
    extendCrawl();
}
