/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef CXSTORAGEDEVICE_H
#define CXSTORAGEDEVICE_H

#include "CxTypes.h"
#include <map>
#include <vector>

namespace ct {

typedef uint64_t FOffset;
uint64_t const UnknownSize = static_cast<FOffset>(-1);

struct FStorageBlock;
struct FMemoryStack;

/// An object managing persistent storage. Abstracts real or virtual file
/// systems. Intended for compatibility with record-oriented integral
/// interfaces.
struct FStorageDevice {
  // identifies storage records within the device.
  // Resource is owned by the device, not by the FRecord structure.
  struct FRecord {
    uint iFile, iRecord;
    FOffset BaseOffset, EndOffset;
    uint Id; // for internal use of storage device
    FRecord() : BaseOffset(0), EndOffset(UnknownSize), Id(0) {}
    FRecord(uint iFile_, uint iRecord_ = 0, FOffset nBaseOffset_ = 0,
            FOffset nLength_ = UnknownSize)
        : iFile(iFile_), iRecord(iRecord_), BaseOffset(nBaseOffset_),
          EndOffset(UnknownSize), Id(0) {
      if (nLength_ != UnknownSize)
        SetLength(nLength_);
    }

    void SetLength(FOffset NewLength) const;
  };

  // allocate a new record for temporary data.
  virtual FRecord AllocNewRecord(FOffset SizeInBytes) = 0;
  FStorageBlock AllocNewBlock(FOffset SizeInBytes);

  virtual void Write(FRecord const &r, void const *pData, FOffset nLength,
                     FOffset Offset = 0) const = 0;
  virtual void Read(FRecord const &r, void *pData, FOffset nLength,
                    FOffset Offset = 0) const = 0;
  virtual void Reserve(FRecord const &r, FOffset nLength) const = 0;
  virtual void Delete(FRecord const &r) = 0;

  virtual ~FStorageDevice(); // = 0;

  FStorageDevice(){};

private:
  FStorageDevice(FStorageDevice const &other); // not implemented, don't copy!
  void operator=(FStorageDevice const &other); // not implemented, don't copy!
};

typedef FStorageDevice::FRecord FRecord;

/// Convenience-combination of storage device and block. Note that these
/// objects DO NOT OWN the data they refer to. The data belongs to the storage
/// device. As such storage blocks are copy-able and don't contain state.
struct FStorageBlock {
  FStorageDevice *pDevice;
  FStorageDevice::FRecord Record;

  FStorageBlock() : pDevice(0){};
  FStorageBlock(FStorageDevice &Device_, FStorageDevice::FRecord const &Record_)
      : pDevice(&Device_), Record(Record_){};

  void RawWrite(void const *pData, FOffset nLength, FOffset Offset = 0) {
    pDevice->Write(Record, pData, nLength, Offset);
  }
  void RawRead(void *pData, FOffset nLength, FOffset Offset = 0) const {
    pDevice->Read(Record, pData, nLength, Offset);
  }

  template <class FScalar>
  void Write(FScalar const *pData, FOffset nLength, FOffset Offset = 0) {
    RawWrite(pData, sizeof(FScalar) * nLength, sizeof(FScalar) * Offset);
  }

  template <class FScalar>
  void Read(FScalar *pData, FOffset nLength, FOffset Offset = 0) const {
    RawRead(pData, sizeof(FScalar) * nLength, sizeof(FScalar) * Offset);
  }

  void Delete() {
    if (pDevice != 0)
      pDevice->Delete(Record);
  }
};

std::ostream &operator<<(std::ostream &out, FStorageDevice::FRecord const &r);
std::ostream &operator<<(std::ostream &out, FStorageBlock const &r);

/// Implements FStorageDevice interface based on basic posix file system
/// interfaces. Each record is represented by a file; files are put into
/// locations returned by ::tmpfile().
struct FStorageDevicePosixFs : public FStorageDevice {
  FRecord AllocNewRecord(FOffset SizeInBytes); // override.
  void Write(FRecord const &r, void const *pData, FOffset nLength,
             FOffset Offset = 0) const; // override
  void Read(FRecord const &r, void *pData, FOffset nLength,
            FOffset Offset = 0) const;                   // override
  void Reserve(FRecord const &r, FOffset nLength) const; // override
  void Delete(FRecord const &r);                         // override

  FStorageDevicePosixFs();
  ~FStorageDevicePosixFs();

protected:
  typedef std::map<uint, FILE *> FFileIdMap;
  FFileIdMap FileIds;
  FILE *GetHandle(FRecord const &r) const;
};

/// Implements FStorageDevice with data kept in memory.
/// Each record is represented by a in-memory storage block.
struct FStorageDeviceMemoryBuf : public FStorageDevice {
  FRecord AllocNewRecord(FOffset SizeInBytes = 0);
  void Write(FRecord const &r, void const *pData, FOffset nLength,
             FOffset Offset = 0) const;
  void Read(FRecord const &r, void *pData, FOffset nLength,
            FOffset Offset = 0) const;
  void Reserve(FRecord const &r, FOffset nLength) const; // override
  void Delete(FRecord const &r) const;

  virtual ~FStorageDeviceMemoryBuf();

private:
  struct FBuffer {
    char *p;
    std::size_t Length;
  };
  typedef std::vector<FBuffer> FBufferSet;
  FBufferSet m_Buffers;
  FBuffer &GetBuf(FRecord const &r) const;
};

enum FFileAndMemOp { OP_AddToMem, OP_Dot, OP_AddToFile, OP_WriteFile };

// provides some simple constant buffer size file+memory operations (see .cpp).
// This template is bound for 64bit and 32bit floats.
template <class FScalar>
void FileAndMemOp(FFileAndMemOp Op, FScalar &f, FScalar *pMemBuf,
                  std::size_t nMemLength, FStorageBlock const &Rec,
                  FMemoryStack &Temp);

} // namespace ct

#endif // CXSTORAGEDEVICE_H
