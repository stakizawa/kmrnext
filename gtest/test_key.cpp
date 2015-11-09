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
      array_key1_ = new size_t[key_size_];
      for (size_t i = 0; i < key_size_; i++) {
	array_key1_[i] = key_size_ - (i + 1);
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
    size_t *array_key0_;    // {0,1,2}
    size_t *array_key1_;    // {2,1,0}
  };

  TEST_F(KeyTest, Constructor) {
    kmrnext::Key k(key_size_);
    EXPECT_EQ(key_size_, k.size());
    EXPECT_THROW({kmrnext::Key ke(key_size_error_);}, std::runtime_error);
  }

  TEST_F(KeyTest, Dim) {
    kmrnext::Key k(key_size_);
    k.set(array_key0_);
    EXPECT_EQ(array_key0_[0], k.dim(0));
    EXPECT_EQ(array_key0_[1], k.dim(1));
    EXPECT_EQ(array_key0_[2], k.dim(2));
    EXPECT_THROW({k.dim(3);}, std::runtime_error);
    EXPECT_THROW({k.dim(key_size_error_);}, std::runtime_error);
  }

  TEST_F(KeyTest, Set_dim) {
    kmrnext::Key k(key_size_);
    k.set(array_key0_);
    k.set_dim(0, 100);
    EXPECT_EQ((size_t)100, k.dim(0));
    EXPECT_EQ(array_key0_[1], k.dim(1));
    EXPECT_EQ(array_key0_[2], k.dim(2));
    EXPECT_THROW({k.set_dim(3, 0);}, std::runtime_error);
    EXPECT_THROW({k.set_dim(key_size_error_, 0);}, std::runtime_error);
  }

  TEST_F(KeyTest, To_string) {
    kmrnext::Key k(key_size_);
    k.set(array_key0_);
    std::string kstr = k.to_string();
    EXPECT_STREQ("<0,1,2>", kstr.c_str());
  }

  TEST_F(KeyTest, Equal) {
    kmrnext::Key k0(key_size_);
    kmrnext::Key k1(key_size_);
    kmrnext::Key k2(key_size_);
    k0.set(array_key0_);
    k1.set(array_key0_);
    k2.set(array_key1_);
    EXPECT_EQ(k1, k0);
    EXPECT_NE(k2, k0);
  }
}  // namespace
