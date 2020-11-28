#include "pch.h"
#include "../GeneralFunctions.h"

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
