#ifndef EXACTLCPK_H
#define EXACTLCPK_H

#include "defs.hpp"
#include "AppConfig.hpp"
#include "ReadsDB.hpp"
#include "rmq_support_sparse_table.hpp"

struct InternalNode{
    int32_t m_leftBound;
    int32_t m_rightBound;
    int32_t m_stringDepth;
    int32_t m_delta;

    bool operator< (const InternalNode& other) const {
        return
            (m_leftBound < other.m_leftBound) ? true :
            ((m_leftBound == other.m_leftBound) ? (m_rightBound < other.m_rightBound) :
             false);
    }

    bool operator== (const InternalNode& other) const{
        return
            (m_leftBound == other.m_leftBound) &&
            (m_rightBound == other.m_rightBound) &&
            (m_stringDepth == other.m_stringDepth);
    }

    InternalNode(){
        m_leftBound = m_rightBound = -1;
        m_stringDepth = m_delta = 0;
    }

    void write(std::ostream& ots, const char *sepStr = "\t") const{
        ots << m_leftBound << sepStr << m_rightBound << sepStr
            << m_stringDepth << sepStr << m_delta;
    };

    void writeln(std::ostream& ots) const{
        write(ots);
        ots << std::endl;
    }

    void dwriteln(std::ostream& ots) const{
        ots << " [";
        write(ots, ",\t");
        ots << "]," << std::endl;
    }

};

struct L1Suffix{
    int32_t m_startPos; // starting position
    int32_t m_errSAPos; // position after one error's SA loc.
    int32_t m_srcStr;   // source string

    bool operator< (const L1Suffix& other) const {
        return( m_errSAPos < other.m_errSAPos );
    }

    L1Suffix(int32_t spos, int32_t epos, int32_t src):
        m_startPos(spos), m_errSAPos(epos), m_srcStr(src) { }

    L1Suffix(){}

    void write(std::ostream& ots, const char *sepStr = "\t") const{
        ots << m_startPos << sepStr << m_errSAPos << sepStr
            << m_srcStr;
    }

    void emit(std::ostream& ots) const{
        ots << m_startPos << "\t" << m_srcStr;
    }

    void writeln(std::ostream& ots) const{
        write(ots);
        ots << std::endl;
    }

    void dwrite(std::ostream& ots) const{
        ots << "\t[";
        write(ots, ",\t");
        ots << "],";
    }

    void dwriteln(std::ostream& ots) const{
        dwrite(ots);
        ots << std::endl;
    }

};

class UpperBoundCheck{
public:
    bool operator()(const int32_t& value, const int32_t& bound) const{
        return (value < bound);
    }
};

class LowerBoundCheck{
public:
    bool operator()(const int32_t& value, const int32_t& bound) const{
        return (value > bound);
    }
};

class IncrPointer{
public:
    void operator()(int32_t& value){
        value += 1;
    }
};

class DecrPointer{
public:
    void operator()(int32_t& value){
        value -= 1;
    }
};

class ExactLCPk{
private:
    AppConfig& m_aCfg;
    int32_t m_strLengths[2];
    int32_t m_strLenPfx[2];
    int32_t m_shiftPos[2];
    std::string m_strXY;
    ivec_t m_gsa, m_gisa, m_glcp;
    ivec_t m_klcpXY[2][2];
    int m_kv;
    double m_nPass;
    double m_passSizes;

    rmq_support_sparse_table<ivec_t, true, ivec_t> m_rangeMinQuery;

    int32_t leftBound0(int32_t curLeaf);
    int32_t rightBound0(int32_t curLeaf);
    void selectInternalNodes0(std::vector<InternalNode>& uNodes);
    void chopPrefix0(const InternalNode& uNode,
                     std::vector<L1Suffix>& leaves);

    void updateExactLCPk(InternalNode& uNode, std::vector<L1Suffix>& leaves);

    void eliminateDupes(std::vector<InternalNode>& uNodes);

    inline int32_t rangeMinLCP(const int32_t& t1, const int32_t& t2){
        if(t1 < 0 || t2 < 0)
            return 0;
        int32_t mxv = std::max(t1, t2);
        if(mxv >= (int32_t)m_gsa.size())
            return 0;
        int32_t mnv = std::min(t1, t2);
        int32_t rpos = m_rangeMinQuery(mnv == mxv ? mnv : (mnv + 1), mxv);
        return m_glcp[rpos];
    }

    inline int32_t rangeMinLCP(const L1Suffix& m1, const L1Suffix& m2){
        return rangeMinLCP(m1.m_errSAPos, m2.m_errSAPos);
    }

    inline int32_t sufRangeMinLCP(const int32_t& t1, const int32_t& t2){
        if(t1 < 0 || t2 < 0)
            return 0;
        if(std::max(t1, t2) >= (int32_t)m_gsa.size())
            return 0;
        int32_t st1 = m_gisa[t1];
        int32_t st2 = m_gisa[t2];
        int32_t mxv = std::max(st1, st2);
        int32_t mnv = std::min(st1, st2);
        int32_t rpos = m_rangeMinQuery(mnv == mxv ? mnv : (mnv + 1), mxv);
        return m_glcp[rpos];
    }

    inline int32_t updatePassLCP(const int32_t& t1, const int32_t& t2){
        if(t1 < 0 || t2 < 0)
            return 0;
        int32_t mxv = std::max(t1, t2);
        if(mxv > (int32_t)m_gsa.size())
            return 0;
        int32_t mnv = std::min(t1, t2);
        int32_t rpos = m_rangeMinQuery(mnv == mxv ? mnv : (mnv + 1), mxv);
        return 1 + m_glcp[rpos];
    }

    inline int32_t updatePassLCP(const L1Suffix& m1, const L1Suffix& m2){
        return updatePassLCP(m1.m_errSAPos, m2.m_errSAPos);
    }

    inline int32_t strPos(const InternalNode& uNode, const L1Suffix& sfx){
        return sfx.m_startPos - uNode.m_delta - m_shiftPos[sfx.m_srcStr];
    }

    template<typename BoundChecker, typename NextPointer>
    void updatePass(int32_t src_ptr, int32_t tgt_ptr,
                    const int32_t& tgt_bound,
                    const InternalNode& uNode,
                    const std::vector<L1Suffix>& leaves
#ifdef DEBUG
                    ,
                    const std::string& dbgStr
#endif
                    ) {
        BoundChecker bound_check;
        NextPointer next_ptr;
        // move the pointer until we reach the first src, target
        while(bound_check(tgt_ptr, tgt_bound)){
            if(leaves[tgt_ptr].m_srcStr != leaves[src_ptr].m_srcStr)
                break;
            src_ptr = tgt_ptr;
            next_ptr(tgt_ptr);
        }
        // if not within the bounds leave
        if(!bound_check(tgt_ptr, tgt_bound) ||
           leaves[tgt_ptr].m_srcStr == leaves[src_ptr].m_srcStr)
            return;
        while(true){
            int32_t rmin = 0;
            int32_t tgt = 0,
                tpos = strPos(uNode, leaves[tgt_ptr]);
            // - get LCP between src_ptr and tgt_ptr from RMQ
                rmin = updatePassLCP(leaves[src_ptr], leaves[tgt_ptr]);
            int32_t score = uNode.m_stringDepth + uNode.m_delta + rmin;
#ifdef DEBUG
            m_aCfg.lfs << "\t[ \""
                       << dbgStr << "\",\t"
                       << score << ",\t"
                       << tgt << ",\t"
                       << tpos << ",\t"
                       << uNode.m_stringDepth << ",\t"
                       << uNode.m_delta << ",\t"
                       << rmin << ",\t"
                ;
            leaves[src_ptr].write(m_aCfg.lfs, ",\t");
            m_aCfg.lfs << ",\t";
            leaves[tgt_ptr].write(m_aCfg.lfs, ",\t");
            m_aCfg.lfs  << "]," << std::endl;
#endif
            assert(tpos >= 0);
            assert(tpos < (int32_t)m_klcpXY[0][1].size());
            // - update target's LCP, if score is higher
            if(score > m_klcpXY[0][1][tpos]){
                int32_t pos = strPos(uNode, leaves[src_ptr]);
                if(pos != tpos){
                    m_klcpXY[0][0][tpos] = pos;
                    m_klcpXY[0][1][tpos] = score;
                }
            }
            // - update tgt_ptr; quit if out of bounds
            int32_t prev_tgt = tgt_ptr;
            next_ptr(tgt_ptr);
            if(!bound_check(tgt_ptr, tgt_bound))
                break;
            // - update src_ptr, if tgt_ptr switches string source
            if(leaves[tgt_ptr].m_srcStr == leaves[src_ptr].m_srcStr)
                src_ptr = prev_tgt;
        }
    }

    void computeK(InternalNode& uNode, std::vector<L1Suffix>& uLeaves,
                  int searchLevel);
    void computeK();
    int32_t leftBoundK(const std::vector<L1Suffix>& trieLeaves,
                       int32_t curLeaf);
    int32_t rightBoundK(const std::vector<L1Suffix>& trieLeaves,
                        int32_t curLeaf);
    void selectInternalNodesK(const InternalNode& prevNode,
                              const std::vector<L1Suffix>& leaves,
                              std::vector<InternalNode>& trieNodes);
    void chopPrefixK(const InternalNode& uNode,
                     const std::vector<L1Suffix>& uLeaves,
                     std::vector<L1Suffix>& trieLeaves);
    void compute0();
    void selectSuffixes0(const InternalNode& uNode,
                         std::vector<L1Suffix>& leaves);
public:
    friend class HeuristicLCPk;
    ExactLCPk(const std::string& x, const std::string& y,
              AppConfig& cfg);
    void print(std::ostream& ofs);
    void compute();
    auto getkLCP() -> const ivec_t (&)[2][2] {
        return m_klcpXY;
    }
    void computeTest(int k);
};


#endif /* EXACTLCPK_H */
