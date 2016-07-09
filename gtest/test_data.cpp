// This test tests Data class.
#include <gtest/gtest.h>
#include "kmrnext.hpp"

namespace {

  class DataTest : public ::testing::Test {
  protected:
    DataTest() {
      value_ = 1;
      value_size_ = sizeof(int);
      data_ = new kmrnext::Data(&value_, value_size_);
    }

    virtual ~DataTest() {
      delete data_;
    }

    virtual void SetUp() {
      // Code here will be called immediately after the constructor (right
      // before each test).
    }

    virtual void TearDown() {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }

    int value_;            // 1
    size_t value_size_;    // sizeof(int)
    kmrnext::Data *data_;  // Data(1, value_size_)
  };

  TEST_F(DataTest, Constructor) {
    kmrnext::Data d(&value_, value_size_);
    EXPECT_EQ(&value_, d.value());
    EXPECT_EQ(value_size_, d.size());
  }

  TEST_F(DataTest, Set_value) {
    kmrnext::Data d0;
    d0.set_value(*data_);
    // value is same
    EXPECT_EQ(*(int*)data_->value(), *(int*)d0.value());
    // address is not same
    EXPECT_NE(data_->value(), d0.value());
    EXPECT_EQ(data_->size(), d0.size());

    // If *value is already set, it throws a runtime_error.
    kmrnext::Data d1(&value_, value_size_);
    EXPECT_THROW({d1.set_value(*data_);}, std::runtime_error);
  }
}
