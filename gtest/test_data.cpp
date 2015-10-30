// This test tests Data class.
#include <gtest/gtest.h>
#include "kmrnext.hpp"

namespace {

  class DataTest : public ::testing::Test {
  protected:
    DataTest() {
      value_ = 1;
      value_size_ = sizeof(int);
      data_ = new Next::Data(&value_, value_size_);
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

    int value_;         // 1
    size_t value_size_; // sizeof(int)
    Next::Data *data_;  // Data(1, value_size_)
  };

  TEST_F(DataTest, Constructor) {
    Next::Data d(&value_, value_size_);
    EXPECT_EQ(&value_, d.value());
    EXPECT_EQ(value_size_, d.size());
  }

  TEST_F(DataTest, Copy_deep) {
    Next::Data d0(NULL, 0);
    d0.copy_deep(*data_);
    // value is same
    EXPECT_EQ(*(int*)data_->value(), *(int*)d0.value());
    // address is not same
    EXPECT_NE(data_->value(), d0.value());
    EXPECT_EQ(data_->size(), d0.size());

    // If *value is already set, it throws a runtime_error.
    Next::Data d1(&value_, value_size_);
    EXPECT_THROW({d1.copy_deep(*data_);}, std::runtime_error);
  }
}
