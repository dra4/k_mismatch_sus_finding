#include "ReadsDB.hpp"
#include "AppConfig.hpp"
#include "ExactLCPk.hpp"
#include "NaiveLCPk.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

class DNA15LCPTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        std::vector<std::string> bfiles;
        bfiles.push_back("t1.fa");
        bfiles.push_back("t2.fa");
        pacfg1 = new AppConfig(bfiles, "psa.out", "psa.log", 1, 0);
        pacfg2 = new AppConfig(bfiles, "pnaive.out", "pnaive.log", 1, 1);
        prdb = new ReadsDB(pacfg1->dir, pacfg1->ifiles, 0);
        const std::string& sx = prdb->getReadById(0);
        const std::string& sy = prdb->getReadById(1);

        pelxy = new ExactLCPk(sx, sy, *pacfg1);
        pnlxy = new NaiveLCPk(sx, sy, *pacfg2);
    }

    virtual void TearDown() {
        delete pelxy;
        delete pnlxy;
        delete pacfg1;
        delete pacfg2;
        delete prdb;
    }

    ExactLCPk *pelxy;
    NaiveLCPk *pnlxy;
    AppConfig *pacfg1, *pacfg2;
    ReadsDB *prdb;

};

TEST_F(DNA15LCPTest, DNA15LCPZeroCheck){
    pelxy->computeTest(0);
    pnlxy->computeTest(0);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(DNA15LCPTest, DNA15LCPOneCheck){
    pelxy->computeTest(1);
    pnlxy->computeTest(1);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(DNA15LCPTest, DNA15LCPTwoCheck){
    pelxy->computeTest(2);
    pnlxy->computeTest(2);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}

TEST_F(DNA15LCPTest, DNA15LCPThreeCheck){
    pelxy->computeTest(3);
    pnlxy->computeTest(3);
    ASSERT_EQ(pelxy->getkLCP()[0][0].size(), pnlxy->getkLCP()[0][0].size());
    ASSERT_EQ(pelxy->getkLCP()[1][0].size(), pnlxy->getkLCP()[1][0].size());
    ASSERT_EQ(pelxy->getkLCP()[0][1], pnlxy->getkLCP()[0][1]);
    ASSERT_EQ(pelxy->getkLCP()[1][1], pnlxy->getkLCP()[1][1]);
}
