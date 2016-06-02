/** 
 * \file      file.h
 * \brief     File hanlder
 * \details   This header includes entropy function.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __FILEIO_H__
#define __FILEIO_H__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>

/**
 * File IO handler.
 */
namespace io
{
	inline std::string getFileExtension( std::string filename ) 
	{
		return filename.substr( filename.find_last_of(".") + 1 );
	}

/* 	/\** Check file suffix to determin gz file. *\/ */
/* 	bool isGzipFile( const char *filename ) */
/* 	{ */
/* 		std::string str = std::string(filename); */
/* 		return isGzipFile( str ); */
/* 	} */

/* 	bool isGzipFile( std::string &filename ) */
/* 	{ */
/* 		if ( filename.substr( filename.find_last_of(".") + 1 ) == "gz" ) */
/* 			return true; */
/* 		return false; */
/* 	} */

	/** Open file. */
	inline void openFile( std::fstream &file, 
						  const char* name,
						  std::_Ios_Openmode mode )
	{
		file.open(name, mode);
		if ( !file) {
			std::cerr << "Can't open " << name << "\n";
			exit (1);
		}
	}
//
//	/** Extract gz format file. */
//	inline void expandFile(const char* input, const char* output)
//	{
//		gzFile infile = gzopen(input, "rb");
//		FILE *outfile = fopen(output, "wb");
//
//		if (!infile ||!outfile) exit(1);
//
//		char buffer[128];
//		int num_read = 0;
//		while ((num_read = gzread(infile, buffer, sizeof(buffer))) > 0) {
//			fwrite(buffer, 1, num_read, outfile);
//		}
//		gzclose(infile);
//		fclose(outfile);
//	}
//
//	/** Is gz file? */
//	//-------------------------------------------------
//	// Very naive file checking only based on file name.
//	// Do robustic thing later
//	//-------------------------------------------------
//	inline bool isCompressed( std::string input_file )
//	{
//			std::string ext, infile, tmpfile;
//		int indexs, indexd;
//
//		infile = input_file;
//		indexd = infile.find_last_of(".");
//		indexs = infile.find_last_of("/");
//
//		ext = infile.substr(indexd+1);
//
//		if ( ext == "gz" || ext == "gunzip" ) {
//			return true;
//		}
//		return false;
//	}
//
//	/** Uncompress a file */
//	inline std::string uncompress( std::string input_file,
//							std::string &tmpdir )
//	{
//		std::string ext, infile, tmpfile;
//		int indexs, indexd;
//
//		infile = input_file;
//		indexd = infile.find_last_of(".");
//		indexs = infile.find_last_of("/");
//
//		ext = infile.substr(indexd+1);
//
//		tmpfile = tmpdir + "/" + infile.substr( indexs+1, indexd-indexs-1 );
//		io::expandFile( infile.c_str(), tmpfile.c_str() );
//
//		return tmpfile;
//	}

//	/** Remove file */
//	inline void removeFile(const char* file)
//	{
//		if( remove(file) != 0 )
//			perror( "Error deleting file" );
//	}

//	/** Umcompress files */
//	inline std::vector<std::string> uncompressFiles( std::vector<std::string> &input_files,
//													 std::string &tmp_dir )
//	{
//		size_t nfile = input_files.size();
//		std::vector<std::string> expanded(nfile);
//		for ( size_t i = 0; i < nfile; i++ )
//			expanded[i] = io::uncompress(input_files[i], tmp_dir);
//
//		return expanded;
//	}
//
//	/** Remove files */
//	inline void removeFiles(std::vector<std::string>& input_files)
//	{
//		for ( size_t i = 0; i < input_files.size(); i++ )
//			io::removeFile(input_files[i].c_str());
//	}
}

#endif
