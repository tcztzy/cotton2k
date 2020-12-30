#include "pch.h"
#include "GeneralFunctions.h"

TEST(GeneralFunctions, LeapYear) {
    EXPECT_EQ(LeapYear(1900), 0);
    EXPECT_EQ(LeapYear(2000), 1);
    EXPECT_EQ(LeapYear(2019), 0);
    EXPECT_EQ(LeapYear(2020), 1);
    EXPECT_EQ(LeapYear(2021), 0);
}

TEST(GeneralFunctions, DateToDoy) {
    EXPECT_EQ(DateToDoy("01-OCT-2020", 2020), 275);
}

TEST(GeneralFunctions, DoyToDate) {
    EXPECT_EQ(DoyToDate(275, 2020), "01-OCT-2020");
}
double SoilTemOnRootGrowth(double);
TEST(RootGrowth, SoilTemOnRootGrowth) {
    EXPECT_EQ(SoilTemOnRootGrowth(30), 1);
    EXPECT_NEAR(SoilTemOnRootGrowth(28), 0.9712, 1e-5);
    EXPECT_NEAR(SoilTemOnRootGrowth(26), 0.9168, 1e-5);
    EXPECT_NEAR(SoilTemOnRootGrowth(24), 0.8368, 1e-5);
    EXPECT_NEAR(SoilTemOnRootGrowth(22), 0.7312, 1e-5);
    EXPECT_NEAR(SoilTemOnRootGrowth(20), 0.6000, 1e-5);
    EXPECT_NEAR(SoilTemOnRootGrowth(18), 0.4432, 1e-5);
    EXPECT_NEAR(SoilTemOnRootGrowth(16), 0.2608, 1e-5);
    EXPECT_NEAR(SoilTemOnRootGrowth(14), 0.0528, 1e-5);
    EXPECT_EQ(SoilTemOnRootGrowth(13.5), 0);
}
