#include "NaiveLCPk.hpp"

NaiveLCPk::NaiveLCPk(const std::string& x, const std::string& y,
                     AppConfig& cfg): sx(x), sy(y), m_aCfg(cfg){
    m_kv = m_aCfg.kv;
}

void NaiveLCPk::runLCPk(const std::string& sa, const std::string& sb,
                        int tidx){
    for(unsigned astart = 0; astart < sa.size(); astart++){
        for(unsigned bstart = 0; bstart < sb.size(); bstart++){
            unsigned howfar = 0, kdx = 0;
            while((bstart + howfar < sb.size()) &&
                  (astart + howfar < sa.size())){
                if(sa[astart + howfar] != sb[bstart + howfar])
                    kdx += 1;
                if(m_kv < (int)kdx)
                    break;
                howfar++;
            }
            if(m_klcpXY[tidx][1][astart] < (int32_t)howfar){
                m_klcpXY[tidx][0][astart] = (int32_t)bstart;
                m_klcpXY[tidx][1][astart] = (int32_t)howfar;
            }
        }
    }
}

void NaiveLCPk::compute(){
    int32_t strLengths[2];
    strLengths[0] = sx.size(); strLengths[1] = sy.size();

    // resize the LCP arrays
    for(unsigned i = 0; i < 2; i++){
        for(unsigned j = 0; j  < 2;j++){
            m_klcpXY[i][j].resize(strLengths[i], 0);
        }
    }
    runLCPk(sx, sy, 0);
    runLCPk(sy, sx, 1);
}

void NaiveLCPk::computeTest(int k){
    m_kv = k;
    compute();
}
