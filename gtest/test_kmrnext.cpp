// This test tests KMRNext class.
#include <gtest/gtest.h>
#include "kmrnext.hpp"

extern kmrnext::KMRNext *gNext;

namespace {

  // The fixture for testing class KMRNext.
  class KMRNextTest : public ::testing::Test {
    virtual void SetUp() {
      gNext->set_io_mode(kmrnext::Memory);
    }

    virtual void TearDown() {
      gNext->set_io_mode(kmrnext::Memory);
    }
  };

  TEST_F(KMRNextTest, IO_mode) {
    EXPECT_EQ(kmrnext::Memory, gNext->io_mode());
    kmrnext::DataStore *ds0 = gNext->create_ds(2);
    EXPECT_EQ(kmrnext::Memory, ds0->io_mode());
    delete ds0;

    gNext->set_io_mode(kmrnext::File);
    EXPECT_EQ(kmrnext::File, gNext->io_mode());
    kmrnext::DataStore *ds1 = gNext->create_ds(2);
    EXPECT_EQ(kmrnext::File, ds1->io_mode());
    delete ds1;
  }
}
