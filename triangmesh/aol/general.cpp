#include <general.h>

namespace aol {

string strprintf(const char * format, ...) {
  // declare variable argument list
  va_list az;
  // copy from my input into this list(second argument is nothing really used, but the last named argument of myself)
  va_start(az, format);
  // give this argument list variable to vscprintf instead of my own arguments (that is in what functions like vsprintf differ from functions like sprintf)
  const int sizeNeeded = vsnprintf ( NULL, 0, format, az ) + 1;
  // restore stack into clean state:
  va_end(az);

  char *buffer = new char[sizeNeeded];

  va_start(az, format);
  vsprintf (buffer, format, az);
  va_end(az);

  // automatic return type conversion into string:
  string ret = buffer;
  delete[] buffer;
  return ret;
}

void makeDirectory ( const char *DirectoryName, bool verbose ) {
  struct stat directory;
  const int statReturn = stat ( DirectoryName, &directory );
  if ( ( statReturn == - 1 ) || !S_ISDIR ( directory.st_mode ) ) {
    string systemCommand = "mkdir \"";
    systemCommand += DirectoryName;
    systemCommand += "\"";
    if ( system ( systemCommand.c_str() ) != EXIT_SUCCESS )
      cerr << "aol::makeDirectory: Calling '" << systemCommand << "' returned an error.\n";
    if ( verbose )
      cerr << "Created directory " << DirectoryName << endl;
  } else {
    if ( verbose )
      cerr << "Directory " << DirectoryName << " already exists\n";
  }
}

bool fileExists ( string filename ) {
  struct stat buf;
  return !stat ( filename.c_str (), &buf ) && S_ISREG ( buf.st_mode );
}

} // namespace aol
