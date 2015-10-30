// This test tests Key and View class.
#include <gtest/gtest.h>
#include "kmrnext.hpp"

namespace {

  // The fixture for testing class Key.
  class KeyTest : public ::testing::Test {
  protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    KeyTest() {
      // You can do set-up work for each test here.
      key_size_ = 3;
      key_size_error_ = 1000;
      array_key0_ = new size_t[key_size_];
      for (size_t i = 0; i < key_size_; i++) {
	array_key0_[i] = i;
      }
    }

    virtual ~KeyTest() {
      // You can do clean-up work that doesn't throw exceptions here.
      delete array_key0_;
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp() {
      // Code here will be called immediately after the constructor (right
      // before each test).
    }

    virtual void TearDown() {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }

    // Objects declared here can be used by all tests in the test case for Key.
    size_t key_size_;       // 3
    size_t key_size_error_; // 1000
    size_t *array_key0_;          // {0,1,2}
  };

  TEST_F(KeyTest, Constructor) {
    Next::Key k(key_size_);
    EXPECT_EQ(key_size_, k.size());
    EXPECT_THROW({Next::Key ke(key_size_error_);}, std::runtime_error);
  }

  TEST_F(KeyTest, Dim) {
    Next::Key k(key_size_);
    k.set(array_key0_);
    EXPECT_EQ(array_key0_[0], k.dim(0));
    EXPECT_EQ(array_key0_[1], k.dim(1));
    EXPECT_EQ(array_key0_[2], k.dim(2));
    EXPECT_THROW({k.dim(3);}, std::runtime_error);
    EXPECT_THROW({k.dim(key_size_error_);}, std::runtime_error);
  }

  TEST_F(KeyTest, Set_dim) {
    Next::Key k(key_size_);
    k.set(array_key0_);
    k.set_dim(0, 100);
    EXPECT_EQ((size_t)100, k.dim(0));
    EXPECT_EQ(array_key0_[1], k.dim(1));
    EXPECT_EQ(array_key0_[2], k.dim(2));
    EXPECT_THROW({k.set_dim(3, 0);}, std::runtime_error);
    EXPECT_THROW({k.set_dim(key_size_error_, 0);}, std::runtime_error);
  }

  TEST_F(KeyTest, To_string) {
    Next::Key k(key_size_);
    k.set(array_key0_);
    std::string kstr = k.to_string();
    EXPECT_STREQ("<0,1,2>", kstr.c_str());
  }
}  // namespace
