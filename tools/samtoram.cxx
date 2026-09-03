#include "ramcore/SamToTTree.h"
#include "ttree/RAMRecord.h"
#include <Compression.h>
#include <cstdio>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <input.sam> [output.root] [-compression N]\n", argv[0]);
        printf("  -compression N  ROOT compression code, algorithm*100+level\n");
        printf("                  (207 = LZMA-7, the default; 505 = ZSTD-5; 404 = LZ4-4; 101 = ZLIB-1)\n");
        return 1;
    }

    const char* input = argv[1];
    const char* output = "ramexample.root";
    // Codes are algorithm*100+level, the same convention the RNTuple tools use.
    // The writer previously took a bare algorithm id and pinned the level to 1,
    // so a TTree file could not be built at the level its RNTuple counterpart used.
    int compression = 207; // LZMA level 7

    for (int i = 2; i < argc; i++) {
        const std::string arg = argv[i];
        if (arg == "-compression") {
            if (i + 1 >= argc) {
                printf("samtoram: -compression needs a value\n");
                return 1;
            }
            compression = std::stoi(argv[++i]);
        } else if (!arg.empty() && arg[0] != '-') {
            output = argv[i];
        }
    }

    samtoram(input, output, true, true, true, compression, RAMRecord::kPhred33);
    
    return 0;
}

