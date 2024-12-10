#include "csv_writer.hh"
#include <filesystem>

namespace df::infra {

  // CSVWriter

  // ctors

  // default: use std::cout
  // Note: Create shared_ptr _os from raw pointer to std::cout;
  //       must use custom deleter (2nd constructor argument)
  //       for shared_ptr to avoid double free (we don't own std::cout).
  //       https://stackoverflow.com/a/17220324
  CSVWriter::CSVWriter()
    : _os(&std::cout, [](std::ostream*){}) {
  }

  CSVWriter::CSVWriter(const std::string& aFile, const bool aAppend, const std::string& aSep)
    : _os(std::make_shared<std::ofstream>(aFile,
                                          aAppend ? std::ofstream::app : std::ofstream::trunc)),
      _sep(aSep) {
    std::filesystem::path lFilePath{aFile};
    std::filesystem::path lDir{std::filesystem::absolute(lFilePath).parent_path()};
    // removing the filename does not work if file is simply "foo.csv", because path{"foo.csv"}.remove_filename() == "".
    //std::filesystem::path lDir{lFilePath.remove_filename()};
    if (!std::filesystem::exists(lDir)) {
      throw std::runtime_error("Directory '" + lDir.string() + "' does not exist");
    }
  }

  CSVWriter::CSVWriter(std::ostream& aStream, const std::string& aSep)
    : _os(&aStream, [](std::ostream*){}),  // see CSVWriter()
      _sep(aSep) {
  }

  // public member functions

  CSVWriter&
  CSVWriter::newline() {
    *_os << NEWLINE;
    if (_flushPolicy == FlushPolicy::kField || _flushPolicy == FlushPolicy::kNewline) {
      _os->flush();
    }
    _colIndex = 0;
    return *this;
  }



}
