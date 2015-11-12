#include <gtest/gtest.h>
#include "kmrnext.hpp"

kmrnext::KMRNext *gNext;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  gNext = kmrnext::KMRNext::init(argc, argv);
  int result = RUN_ALL_TESTS();
  kmrnext::KMRNext::finalize();
  return result;
}
