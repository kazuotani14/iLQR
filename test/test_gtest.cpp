#include "gtest/gtest.h"

TEST(GTestTest, DoesItWork) {
  EXPECT_EQ(5.0, 4.0+1);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
