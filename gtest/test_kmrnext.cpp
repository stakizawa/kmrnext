// This test tests KMRNext class.
#include <gtest/gtest.h>
#include "kmrnext.hpp"

extern kmrnext::KMRNext *gNext;

namespace {

  // The fixture for testing class KMRNext.
  class KMRNextTest : public ::testing::Test {
    virtual void SetUp() {
      gNext->set_io_mode(kmrnext::KMRNext::Memory);
    }

    virtual void TearDown() {
      gNext->set_io_mode(kmrnext::KMRNext::Memory);
    }
  };

  TEST_F(KMRNextTest, IO_mode) {
    EXPECT_EQ(kmrnext::KMRNext::Memory, gNext->io_mode());
    gNext->set_io_mode(kmrnext::KMRNext::File);
    EXPECT_EQ(kmrnext::KMRNext::File, gNext->io_mode());
  }
}
