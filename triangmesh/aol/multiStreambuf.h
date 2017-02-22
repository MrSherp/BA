#ifndef __MULTISTREAMBUF_H
#define __MULTISTREAMBUF_H

#include <general.h>

namespace aol {

// \brief Output stream buffer that sends its content to multiple other stream buffers.
class MultiStreambuf : public streambuf {
public:
  MultiStreambuf() { setp(0, 0); }

  size_t addStreambuf(streambuf * buffer) {
    _streambufPtrList.push_back(buffer);
    return _streambufPtrList.size();
  }
  size_t removeStreambuf(streambuf * buffer) {
    _streambufPtrList.remove(buffer);
    return _streambufPtrList.size();
  }

protected:
  typedef std::char_traits<char> traits_type;
  typedef traits_type::int_type  int_type;

  typedef list<streambuf*> StreambufList;
  StreambufList _streambufPtrList;

  int_type overflow( int_type c) {
    if (!traits_type::eq_int_type(c, traits_type::eof())) {
        StreambufList::iterator iter = _streambufPtrList.begin();
        for (; iter != _streambufPtrList.end(); ++iter)
        (*iter)->sputc(c);
    }
    return traits_type::not_eof(c);
  }
  //! synchronizes all buffers. Returns for failure if at least one participation buffer has failed and returned -1. 
  int_type sync() {
    int ret = 0;
    StreambufList::iterator iter = _streambufPtrList.begin();
    for (; iter != _streambufPtrList.end(); ++iter)
        // if ret has already set to value "-1" (failed), do not overwrite this, but keep.
        if ( ret == -1 )
            (*iter)->pubsync();
        // otherwise, give *iter a chance to indicate failure.
        else
            ret = (*iter)->pubsync();
        return ret;
  }
};

//! \brief Print cout and cerr not only to consolen (or whereever you have redirected it), but also to a file.
//  When you create an instance of this class, the standard cout / cerr buffers are replaced by"output buffer collections" that distribute every write attemp to the previously used cout / cerr buffer 
//  as well as to a file output buffer.  It is possible to iterate this procedure by creating multiple instances of AdditionalOutputToFile (but they have to be destroyed in reverse order).
class AdditionalOutputToFile {
public:
  AdditionalOutputToFile ( string Filename )
    : _filename ( Filename )
    , _filestream ( Filename.c_str() )
    , _previousCoutStreambuf ( cout.rdbuf() )
    , _previousCerrStreambuf ( cerr.rdbuf() )
    , _previousClogStreambuf ( clog.rdbuf() ) {

    _multiCoutStreambuf.addStreambuf ( _previousCoutStreambuf ); _multiCerrStreambuf.addStreambuf ( _previousCerrStreambuf ); _multiClogStreambuf.addStreambuf ( _previousClogStreambuf );
    _multiCoutStreambuf.addStreambuf ( _filestream.rdbuf() ); _multiCerrStreambuf.addStreambuf ( _filestream.rdbuf() ); _multiClogStreambuf.addStreambuf ( _filestream.rdbuf() );

    cout.rdbuf ( &_multiCoutStreambuf );
    cerr.rdbuf ( &_multiCerrStreambuf );
    clog.rdbuf ( &_multiClogStreambuf );
  }
  
  ~AdditionalOutputToFile () {
    cout.flush(); cerr.flush(); clog.flush();
    cout.rdbuf ( _previousCoutStreambuf ); cerr.rdbuf ( _previousCerrStreambuf ); clog.rdbuf ( _previousCerrStreambuf );
    clog << "All output has been written to file " << _filename << endl;
  }

protected:
  string      _filename;
  ofstream    _filestream;
  streambuf * _previousCoutStreambuf;
  streambuf * _previousCerrStreambuf;
  streambuf * _previousClogStreambuf;
  MultiStreambuf _multiCoutStreambuf;
  MultiStreambuf _multiCerrStreambuf;
  MultiStreambuf _multiClogStreambuf;
};

} // end of namespace aol.

#endif
