#include "bwa_reader.h"

using namespace std;

BwaReader::BwaReader()
{
    init();
}

BwaReader::BwaReader( string &file, bool p ) 
{
    init();
    bwa_file = file;
	pairend  = p;
}

void BwaReader::init()
{
    verbose     = false;
    set_contig  = false;
    outdir      = ".";
    //aln_offset  = ALLOW;
}

void BwaReader::parse()
{
	double t0 = mytime();
    readFile();
	log("BWA file parsed", t0);
}

void BwaReader::readFile()
{
    string rfile = outdir + "/gabe_reads";
    ofstream out(rfile.c_str(), ios::out);
    if ( !out ) {
        cerr << "Can't open " << rfile << "\n";
        exit (1);
    }

	
	// Skip unmapped reads (0x4) and secondary matches (0x100) and supplementary alignment (0x800)
	bool sam_flag = bwa_file.find( ".sam" ) != std::string::npos ? true : false; 
	std::string cmd = sam_flag ? 
		"samtools view -S -h -F 0x4 " + bwa_file :
		"samtools view    -h -F 0x4 " + bwa_file ;
	if ( pairend ) cmd += " | samtools view -S -h -F 0x8 -";
	cmd += " | samtools view -S -h -F 0x100 - | samtools view -S -h -F 0x800 -";

	std::cerr << "[SAMTOOLS] " << cmd << "\n";

	// Read in streams
	redi::ipstream in( cmd.c_str() );
	
    string line, rname;
    query_t rnum = 1;
	size_t edist;
	out << "#Index\tRead_name\tEdit_distance\n";
	while(getline(in,line)) {
        if ( line[0] == '@' ) getHeader(line);
        else {
            if ( ! set_contig ) setContigs();
            if ( ! parseLine(line, rnum, rname, edist) ) continue;

			out << rnum << "\t" << rname << "\t" << edist << "\n";
			rnum++;			
        }
    }
}

void BwaReader::getHeader( string &line )
{
    istringstream iss(line);
    vector<string> c;
    string tok;
    while(iss>>tok) c.push_back(tok);

    if ( c[0] != "@SQ" ) return;

    // contig name
    size_t pos = c[1].find("SN:");
    assert( pos != string::npos );
    string contig = c[1].substr(3);

    // contig length
    pos = c[2].find("LN:");
    assert( pos != string::npos );
    size_t length = stoul( c[2].substr(3), nullptr, 0 );

    conlens.insert( pair<string, size_t>( contig, length ) );
}

void BwaReader::setContigs()
{
    // map -> ordered
    size_t index = 0;

    for ( auto it = conlens.cbegin(); it != conlens.cend(); ++it, ++index ) {
        string sbjct = it->first;
        size_t length = it->second;
        sbjlens.insert( pair<sbjct_t, size_t>( index, length ) );
        num2con.insert( pair<sbjct_t, string>( index, sbjct ) );
        con2num.insert( pair<string, sbjct_t>( sbjct, index ) );
    }
    conlens.clear();
    set_contig = true;
}

bool BwaReader::parseLine( string &line, query_t rnum, std::string &rname, size_t &edit )
{
    istringstream iss(line);
    vector<string> c;
    string tok;
    while(iss>>tok) c.push_back(tok);

    rname = c[0];
	std::string sbj = c[2];
    if ( sbj == "*" ) return false;

	//------------------------------------------------------------------
    // Check CIGAR string
    //if ( ! goodCigar(c[5]) ) return pair<bool,string>(false, string());
	// Current regex_search is very slow.
	// For now, ignore any clipped, spliced, padded alignments
	//------------------------------------------------------------------
	if ( ! properAlignment(c[5]) ) return false;

	//--------------------------------------------------
    // bwa mem does not have XT: field
    // Therefore, c[12] is not valid for bwa mem - c[11]
	//--------------------------------------------------
    pair<bool,int> srch;
    for ( size_t i = 11; i <= 12; i++ ) {
        //srch = getEditDistance( line );
        srch = getEditDistance( c[i]);
        if ( srch.first ) break;
    }
    if ( !srch.first ) {
        std::cerr << "[Error] edit distance can't be found in column 12 and 13\n"
                  << "Line:" << line << "\n";
        exit(1);
    }

	if ( pairend ) {
		rname += ( stoi(c[1]) & 0x0040 ) ? "/1" : "/2";
	}
	
    edit = srch.second;
	
    SbjctList slist;
    sbjct_t sbjnum = con2num[sbj];
    slist.push_back(sbjnum);
    querys.insert( pair<query_t, SbjctList>( rnum, slist ) );
    auto it = sbjcts.find(sbjnum);
    if ( it == sbjcts.end() )
        sbjcts.insert( pair<sbjct_t, size_t>(sbjnum, 1) );
    else it->second++;
    
    string last = c[c.size()-1];
    if ( last.find( "XA:Z:") == string::npos )
        return true;

    string extra = last.substr(5);
    extra.pop_back(); // remove last ';' character

    iss.clear(); iss.str(extra);
    vector<string> more_hits;
    while(getline(iss,tok,';')) 
        more_hits.push_back(tok);
    
    for ( auto it : more_hits ) {
        istringstream mss(it);
        vector<string> desc;
        while(getline(mss,tok,','))
            desc.push_back(tok);
        
        string more_sbj = desc[0];
        size_t more_edit = stoi(desc.back());
        if ( more_edit > edit ) continue;

        size_t more_num = con2num[more_sbj];
        querys[rnum].push_back( more_num );
        auto mt = sbjcts.find(more_num);
        if ( mt == sbjcts.end() )
            sbjcts.insert( pair<sbjct_t, size_t>(more_num, 1) );
        else
            mt->second++;
    }
	return true;
}

size_t BwaReader::getTypeCount( string &cigar, string &pattern )
{
    size_t count = 0;
    regex reg(pattern);
    sregex_iterator it(cigar.begin(), cigar.end(), reg);
    sregex_iterator end;
    while(it != end) {
        std::smatch m = *it;
        std::string s = m.str();
        s.pop_back();
        size_t n = stoul(s);
        count += n;
        ++it;
    }
    return count;
}

//====================================================================
// regex is slow. So do not use this
//====================================================================
// bool BwaReader::goodCigar( string &cigar )
// {
//     string type = "[SHNP]";
//     regex reg(type);
//     if ( ! regex_search(cigar,reg) ) return true;

//     size_t sum_count = 0;

//     // Soft clipped aligment
//     type = "[0-9]+S";
//     sum_count += getTypeCount( cigar, type );
//     if ( sum_count > aln_offset ) return false;

//     // Hard clipped alignmenbt
//     type = "[0-9]+H";    
//     sum_count += getTypeCount( cigar, type );    
//     if ( sum_count > aln_offset ) return false;
    
//     // Spliced alignment
//     type = "[0-9]+N";    
//     sum_count += getTypeCount( cigar, type );    
//     if ( sum_count > aln_offset ) return false;

//     // Padded alignment
//     type = "[0-9]+P";    
//     sum_count += getTypeCount( cigar, type );    
//     if ( sum_count > aln_offset ) return false;
    
//     return true;
// }


//====================================================================
// Ignore following alignments
// S: soft clipping
// H: hard clipping
// P: padded alignmnet
// N: skpped region from the reference
//====================================================================
bool BwaReader::properAlignment( string &cigar )
{
	for ( size_t i = 0; i < cigar.size(); i++ ) {
		char ch = cigar[i];
		if ( ch == 'S' || 
			 ch == 'H' ||
			 ch == 'N' ||
			 ch == 'P' ) 
			return false;
	}
	return true;
}

pair<bool,int> BwaReader::getEditDistance( string &col )
{
    size_t pos = col.find("NM:i:");
    if ( pos == string::npos ) 
        return pair<bool,int>(false,-1);
        
    assert ( pos != string::npos );
    int dist = stoi(col.substr(5));
    return pair<bool,int>(true,dist);
}

