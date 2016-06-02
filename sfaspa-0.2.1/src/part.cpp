#include "sfa.h"
#include "cmdargs.h"
#include <fstream>
#include <iostream>

size_t getSize( const char *filename, size_t &num_seq )
{
    std::ifstream fstrm(filename, std::ios_base::in | std::ios_base::binary);
    
    if ( !fstrm ) {
        std::cerr << "Can't open " << filename << "\n";
        exit (1);
    }
    
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    if ( fio::getFileExtension( std::string(filename) ) == "gz" ) 
        in.push(boost::iostreams::gzip_decompressor());
    in.push(fstrm);
    std::istream istrm(&in);
    
    size_t len = 0;
    std::string line, tag, seq;
    while ( std::getline(istrm, line) ) {
        if ( line[0] == '>' ) {
            if ( tag != "" && seq != "" ) {
                len += seq.size();
                num_seq++;
            }
            tag = line; seq = ""; 
        }
        else seq += line;
    }

    if ( tag != "" && seq != "") {
        len += seq.size();
        num_seq++;
    }

    fstrm.close();

    return len;
}

int main(int argc, char **argv)
{
    std::vector<std::string> input_files;
    args::parse_cmd_args_parts(argc, argv, input_files);

    args::printCommand(argc, argv, std::cout);

    size_t sum_len = 0;
    size_t max_doc = std::numeric_limits<SfaType>::max();

    size_t num_seq = 0;
    for ( size_t i = 0; i < input_files.size(); i++ ) {
        sum_len += getSize(input_files[i].c_str(), num_seq);
    }
    
    printf("#files:%zu\n", input_files.size());
    printf("#sequences:%zu\n", num_seq);
    printf("#letters:%zu\n", sum_len);
    printf("#maxium:%zu\n", max_doc-num_seq);
    printf("#parts:%d\n", 1 + int(sum_len/(max_doc-num_seq)));

    return 0;
}
