// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause

#ifndef vtkXMLDataHeaderPrivate_DoNotInclude
#error "do not include unless you know what you are doing"
#endif

#ifndef vtkXMLDataHeaderPrivate_h
#define vtkXMLDataHeaderPrivate_h

#include "vtkType.h"
#include <vector>

// Abstract interface using type vtkTypeUInt64 to access an array
// of either vtkTypeUInt32 or vtkTypeUInt64.  Shared by vtkXMLWriter
// and vtkXMLDataParser to write/read binary data headers.

class vtkXMLDataHeader
{
public:
  virtual void Resize(size_t count) = 0;
  virtual vtkTypeUInt64 Get(size_t index) const = 0;
  virtual bool Set(size_t index, vtkTypeUInt64 value) = 0;
  virtual size_t WordSize() const = 0;
  virtual size_t WordCount() const = 0;
  virtual unsigned char* Data() = 0;
  size_t DataSize() const { return this->WordCount() * this->WordSize(); }
  virtual ~vtkXMLDataHeader() = default;
  static inline vtkXMLDataHeader* New(int width, size_t count);
};

template <typename T>
class vtkXMLDataHeaderImpl : public vtkXMLDataHeader
{
  std::vector<T> Header;

public:
  vtkXMLDataHeaderImpl(size_t n)
    : Header(n, 0)
  {
  }
  void Resize(size_t count) override { this->Header.resize(count, 0); }
  vtkTypeUInt64 Get(size_t index) const override { return this->Header[index]; }
  bool Set(size_t index, vtkTypeUInt64 value) override
  {
    this->Header[index] = T(value);
    return vtkTypeUInt64(this->Header[index]) == value;
  }
  size_t WordSize() const override { return sizeof(T); }
  size_t WordCount() const override { return this->Header.size(); }
  unsigned char* Data() override { return reinterpret_cast<unsigned char*>(this->Header.data()); }
};

vtkXMLDataHeader* vtkXMLDataHeader::New(int width, size_t count)
{
  switch (width)
  {
    case 32:
      return new vtkXMLDataHeaderImpl<vtkTypeUInt32>(count);
    case 64:
      return new vtkXMLDataHeaderImpl<vtkTypeUInt64>(count);
  }
  return nullptr;
}

#endif
// VTK-HeaderTest-Exclude: vtkXMLDataHeaderPrivate.h
