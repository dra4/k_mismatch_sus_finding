#include "ReadsDB.hpp"
#include "AppConfig.hpp"
#include "ExactLCPk.hpp"
#include "NaiveLCPk.hpp"
#include "HeuristicLCPk.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

class PrimatesLCPTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        std::vector<std::string> bfiles;
        bfiles.push_back("C.aethiops.fa");
        bfiles.push_back("C.albifrons.fa");
        // bfiles.push_back("t1.fa");
        // bfiles.push_back("t2.fa");
        pacfg1 = new AppConfig(bfiles, "psa.out", "psa.log", 1, 0);
        pacfg2 = new AppConfig(bfiles, "pnaive.out", "pnaive.log", 1, 1);
        pacfg3 = new AppConfig(bfiles, "psa.out", "psa.log", 1, 2, 5);
        pacfg4 = new AppConfig(bfiles, "pnaive.out", "pnaive.log", 1, 2, 5);
        prdb = new ReadsDB(pacfg1->dir, pacfg1->ifiles, 0);
        const std::string& sx = prdb->getReadById(0);
        const std::string& sy = prdb->getReadById(1);

        pelxy = new ExactLCPk(sx, sy, *pacfg1);
        pnlxy = new NaiveLCPk(sx, sy, *pacfg2);
        phlxy1 = new HeuristicLCPk(sx, sy, *pacfg3);
        phlxy2 = new HeuristicLCPk(sx, sy, *pacfg4);
    }

    virtual void TearDown() {
        delete pelxy;
        delete pnlxy;
        delete pacfg1;
        delete pacfg2;
        delete pacfg3;
        delete pacfg4;
        delete phlxy1;
        delete phlxy2;
        delete prdb;
    }

    ExactLCPk *pelxy;
    NaiveLCPk *pnlxy;
    HeuristicLCPk *phlxy1, *phlxy2;
    AppConfig *pacfg1, *pacfg2, *pacfg3, *pacfg4;
    ReadsDB *prdb;

};

TEST_F(PrimatesLCPTest, PrimatesLCPZeroCheck){
    pelxy->computeTest(0);
    pnlxy->computeTest(0);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(PrimatesLCPTest, PrimatesLCPOneCheck){
    pelxy->computeTest(1);
    pnlxy->computeTest(1);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(PrimatesLCPTest, PrimatesLCPTwoCheck){
    pelxy->computeTest(2);
    pnlxy->computeTest(2);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(PrimatesLCPTest, PrimatesLCPThreeCheck){
    pelxy->computeTest(3);
    pnlxy->computeTest(3);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(PrimatesLCPTest, PrimatesHstLCPOne5Check){
    phlxy1->computeTest(1, 5);
    phlxy2->computeCrawlTest(1, 5);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}
