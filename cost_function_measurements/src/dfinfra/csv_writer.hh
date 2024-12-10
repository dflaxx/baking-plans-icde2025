#pragma once

#include "standard_includes.hh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <ostream>


namespace df::infra {

class CSVWriter {
  public:
    enum class FlushPolicy {
      kDefault,  // the ostream's default
      kNewline,  // flush after each newline
      kField     // flush after each field
    };
  public:
    CSVWriter();
    CSVWriter(const std::string& aFile,
                       const bool aAppend = false,
                       const std::string& aSep = ";");
    CSVWriter(std::ostream& aStream,
                       const std::string& aSep = ";");
    CSVWriter(const CSVWriter&) = delete;
    CSVWriter(CSVWriter&&) = delete;

    inline ~CSVWriter() = default;

  public:
    inline void setFlushPolicy(const FlushPolicy aPolicy) { _flushPolicy = aPolicy; }

  public:
    //CSVWriter& writeLine(const std::string& aLine);

    template <typename T>
    inline
    CSVWriter& writeField(const T& aField, const bool aLastInLine = false) {
      if (_colIndex > 0) { *_os << _sep; }
      *_os << aField;
      if (_flushPolicy == FlushPolicy::kField) { _os->flush() ; }
      ++_colIndex;
      if (aLastInLine) { newline(); }
      return *this;
    }

    // with formatting flags from std::ios_base
    template <typename T>
    CSVWriter& writeFieldFmt(const T& aField, const std::ios_base::fmtflags aFlgs);

    // with formatting flags and precision (+ fixed width!)
    template <typename T>
    CSVWriter& writeFieldFmt(const T& aField, const std::ios_base::fmtflags aFlgs, const size_t aPrecision);

    CSVWriter& newline();

  private:
    std::shared_ptr<std::ostream> _os;
    std::string _sep = ";";
    std::string _commentPrefix = "#";  // unused

    const static char NEWLINE = '\n';
    
    uint32_t _colIndex = 0;

    FlushPolicy _flushPolicy = FlushPolicy::kDefault;


};  // CSVWriter

template <typename T>
inline
CSVWriter&
CSVWriter::writeFieldFmt(const T& aField, const std::ios_base::fmtflags aFlgs) {
  // flags vs setf:
  //   std::ios_base::setf(flags): adds flags to existing flags in stream
  //   std::ios_base::flags(flags): overwrites existing flag setting with flags
  const std::ios_base::fmtflags lFlgsBefore = _os->flags();
  _os->setf(aFlgs);
  if (_colIndex > 0) { *_os << _sep; }
  *_os << aField;
  if (_flushPolicy == FlushPolicy::kField) { _os->flush() ; }
  ++_colIndex;
  _os->flags(lFlgsBefore);
  return *this;
}

template <typename T>
inline
CSVWriter&
CSVWriter::writeFieldFmt(const T& aField, const std::ios_base::fmtflags aFlgs, const size_t aPrecision) {
  const std::ios_base::fmtflags lFlgsBefore = _os->flags();
  const size_t lPrecisionBefore = _os->precision();
  _os->setf(aFlgs);
  _os->setf(std::ios_base::fixed);
  _os->precision(aPrecision);
  if (_colIndex > 0) { *_os << _sep; }
  *_os << aField;
  if (_flushPolicy == FlushPolicy::kField) { _os->flush() ; }
  ++_colIndex;
  _os->flags(lFlgsBefore);
  _os->precision(lPrecisionBefore);
  return *this;
}

}  // df::infra
