#ifndef DATA_ELEMENT_HPP
#define DATA_ELEMENT_HPP
/// \file
/// Header that defines interfaces of DataElement classes.

#include "kmrnext.hpp"

namespace kmrnext {

  ///////////////////////////////////////////////////////////////////////////
  /// A class of DataElement that is used in SimpleFileDataStore
  ///////////////////////////////////////////////////////////////////////////
  class SimpleFileDataElement : public DataElement {
    typedef DataElement base;

  public:
    SimpleFileDataElement();

    virtual ~SimpleFileDataElement() {}

    /// It clears all attribute of the Data.
    void clear();

    /// It restores the Data in the DataElement from a file buffer.
    ///
    /// \param[in] buf  A file buffer
    void restore(char* buf);

    /// It is called when the Data in the DataElement is written to a file.
    ///
    /// \param[in] start_pos   File start position
    /// \param[in] written_siz Written size in bytes
    void written(size_t start_pos, size_t written_siz);

    /// It clears memory cache of the DataElement.
    void clear_cache();

  private:
    // True, if Data is updated only on memory
    bool data_updated_;
    // The offset of the Data in a file
    size_t data_file_offset_;
    // The size of the Data in a file
    size_t data_file_size_;

    void set_data(const void* val, const size_t siz, bool overwrite=false);
  };

}

#endif
