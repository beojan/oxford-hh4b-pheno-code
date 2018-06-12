#include <ROOT/TDataFrame.hxx>
#include <ROOT/TCsvDS.hxx>

#include <iostream>
#include <string>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;
namespace TDF = ROOT::Experimental::TDF;
int main(int argc, char* argv[]) {
    if(argc != 2) {
        std::cerr << "Usage: csv2root filename\n";
        return 1;
    }

    auto df = TDF::MakeCsvDataFrame(argv[1]);
    fs::path root_file = argv[1];
    root_file.replace_extension("root");
    df.Snapshot(root_file.stem().native(), root_file.native(), df.GetColumnNames());
    return 0;
}

