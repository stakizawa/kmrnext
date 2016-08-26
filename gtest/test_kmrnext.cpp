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

  TEST_F(KMRNextTest, Init) {
    // Calling KMRNext::init twice or more returns the same KMRNext
    // instance as gNext.
    kmrnext::KMRNext *n0 = kmrnext::KMRNext::init();
    EXPECT_EQ(gNext, n0);
    kmrnext::KMRNext *n1 = kmrnext::KMRNext::init(0, NULL);
    EXPECT_EQ(gNext, n1);
  }
}
