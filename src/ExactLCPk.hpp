#ifndef EXACTLCPK_H
#define EXACTLCPK_H

#include "defs.hpp"
#include "AppConfig.hpp"
#include "ReadsDB.hpp"
#include "rmq_support_sparse_table.hpp"
#include <atomic>

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

    bool operator< (const L1Suffix& other) const {
        return( m_errSAPos < other.m_errSAPos );
    }

    L1Suffix(int32_t spos, int32_t epos):
        m_startPos(spos), m_errSAPos(epos) { }

    L1Suffix(){}

    void write(std::ostream& ots, const char *sepStr = "\t") const{
        ots << m_startPos << sepStr << m_errSAPos;
    }

    void emit(std::ostream& ots) const{
        ots << m_startPos;
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
    int32_t m_strLength;
    std::string m_str;
    ivec_t m_gsa, m_gisa, m_glcp;
    std::atomic<int32_t> *m_klcpXY;
    int m_kv;
    double m_nPass;
    double m_passSizes;
    std::atomic<int32_t> **m_cur_idx_and_limits;

    rmq_support_sparse_table<ivec_t, true, ivec_t> m_rangeMinQuery;

    int32_t leftBound0(int32_t curLeaf);
    int32_t rightBound0(int32_t curLeaf);
    void selectInternalNodes0(std::vector<InternalNode>& uNodes);
    void chopPrefix0(const InternalNode& uNode,
                     std::vector<L1Suffix>& leaves);

    void updateExactLCPk(const InternalNode& uNode, const std::vector<L1Suffix>& leaves);

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
        return sfx.m_startPos - uNode.m_delta;
    }

    void updatePass(const int32_t& tgt_bound,
                    const InternalNode& uNode,
                    const std::vector<L1Suffix>& leaves
#ifdef DEBUG
                    ,
                    const std::string& dbgStr
#endif
                    ) {
        int32_t tgt_ptr = 0;
        const int32_t tgt_bound_minus_one = tgt_bound - 1;
        if(tgt_ptr < tgt_bound_minus_one) {
            int32_t tpos = strPos(uNode, leaves[tgt_ptr]);
            int32_t rmin = updatePassLCP(leaves[tgt_ptr + 1], leaves[tgt_ptr]);
            int32_t score = uNode.m_stringDepth + uNode.m_delta + rmin;
            int32_t observed_prev_score = m_klcpXY[tpos].load();
            while(score > observed_prev_score && !m_klcpXY[tpos].compare_exchange_weak(observed_prev_score, score));
            ++tgt_ptr;
        } else {
            return;
        }
        while(tgt_ptr < tgt_bound_minus_one){
            int32_t tpos = strPos(uNode, leaves[tgt_ptr]);
            int32_t src_ptr = 0;
            int32_t rmin = 0;
            int32_t rmin_a = updatePassLCP(leaves[tgt_ptr - 1], leaves[tgt_ptr]);
            int32_t rmin_b = updatePassLCP(leaves[tgt_ptr + 1], leaves[tgt_ptr]);
            if(rmin_a > rmin_b) {
                src_ptr = tgt_ptr - 1;
                rmin = rmin_a;
            } else {
                src_ptr = tgt_ptr + 1;
                rmin = rmin_b;
            }
            assert(tpos >= 0);
            assert(tpos < (int32_t)m_strLength);
            int32_t score = uNode.m_stringDepth + uNode.m_delta + rmin;
            int32_t observed_prev_score = m_klcpXY[tpos].load();
            while(score > observed_prev_score && !m_klcpXY[tpos].compare_exchange_weak(observed_prev_score, score));
            ++tgt_ptr;
        }
        if(tgt_ptr == tgt_bound_minus_one) {
            int32_t tpos = strPos(uNode, leaves[tgt_ptr]);
            int32_t rmin = updatePassLCP(leaves[tgt_ptr - 1], leaves[tgt_ptr]);
            int32_t score = uNode.m_stringDepth + uNode.m_delta + rmin;
            int32_t observed_prev_score = m_klcpXY[tpos].load();
            while(score > observed_prev_score && !m_klcpXY[tpos].compare_exchange_weak(observed_prev_score, score));
            ++tgt_ptr;
        }
    }

    void launch(const std::vector<InternalNode>& uNodes, const std::vector<int32_t>& indices, size_t tid, unsigned t_count);

    void computeK(const InternalNode& uNode, const std::vector<L1Suffix>& uLeaves,
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
    ExactLCPk(const std::string& s, AppConfig& cfg);
    void print(std::ostream& ofs);
    void compute();
    auto getkLCP() -> std::vector<int32_t> {
        std::vector<int32_t> ret(m_strLength);
        for(size_t i = 0; i < m_strLength; ++i) {
            ret[i] = m_klcpXY[i].load();
        }
        return ret;
    }
    void computeTest(int k);
};


#endif /* EXACTLCPK_H */
