#include "ReadsDB.hpp"
#include "AppConfig.hpp"
#include "HeuristicLCPk.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

class HueristicDNA15LCPTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        std::vector<std::string> bfiles;
        // bfiles.push_back("C.aethiops.fa");
        // bfiles.push_back("C.albifrons.fa");
        bfiles.push_back("t1.fa");
        bfiles.push_back("t2.fa");
        pacfg1 = new AppConfig(bfiles, "psa.out", "psa.log", 1, 2, 5);
        pacfg2 = new AppConfig(bfiles, "pnaive.out", "pnaive.log", 1, 2, 5);
        prdb = new ReadsDB(pacfg1->dir, pacfg1->ifiles, 0);
        const std::string& sx = prdb->getReadById(0);
        const std::string& sy = prdb->getReadById(1);

        pelxy = new HeuristicLCPk(sx, sy, *pacfg1);
        pnlxy = new HeuristicLCPk(sx, sy, *pacfg2);
    }

    virtual void TearDown() {
        delete pelxy;
        delete pnlxy;
        delete pacfg1;
        delete pacfg2;
        delete prdb;
    }

    HeuristicLCPk *pelxy, *pnlxy;
    AppConfig *pacfg1, *pacfg2;
    ReadsDB *prdb;

};

TEST_F(HueristicDNA15LCPTest, HeuristicNaiveLCPZero5Check){
    std::vector<int32_t> xlcp = {8, 7, 6, 9, 3, 2, 8, 8, 3, 4, 5, 4, 3, 2, 1};
    std::vector<int32_t> ylcp = {6, 2, 2, 7, 8, 7, 8, 8, 7, 6, 5, 4, 3, 2, 1};
    pnlxy->computeTest(0, 5);
    ASSERT_EQ(pnlxy->getkLCP()[0][0].size(), xlcp.size());
    ASSERT_EQ(pnlxy->getkLCP()[1][0].size(), ylcp.size());
    ASSERT_EQ(pnlxy->getkLCP()[0][1], xlcp);
    ASSERT_EQ(pnlxy->getkLCP()[1][1], ylcp);
}

TEST_F(HueristicDNA15LCPTest, DNA15LCPZero0Check){
    pnlxy->computeCrawlTest(0, 0);
    pelxy->computeTest(0, 0);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(HueristicDNA15LCPTest, DNA15LCPZero5Check){
    pelxy->computeTest(0, 5);
    pnlxy->computeCrawlTest(0, 5);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(HueristicDNA15LCPTest, DNA15LCPOne5Check){
    pelxy->computeTest(1, 5);
    pnlxy->computeCrawlTest(1, 5);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}
