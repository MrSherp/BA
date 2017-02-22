#include <PParser.h>

namespace aol {

// Class PParser

PParser::PParser () : _values () {}

PParser::PParser ( const PParser& parser ) : _values ( parser._values ) {}

PParser& PParser::operator = ( const PParser& parser ) {
  // Beware of self-assignment
  if ( this == &parser ) return *this;
  _values = parser._values;
  return *this;
}

PParser::~PParser () {}

PParser::PParser ( const std::string& fileName ) : _values () { read ( fileName );}

void PParser::read ( const std::string& fileName ) {
  std::ifstream stream;
  stream.open ( fileName.c_str() );
  if ( stream.fail() ) throw std::invalid_argument ( aol::strprintf ( "cannot open file in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  // Extend this for vadditional types
  std::map <string, TypePtr> types;
  types["int"]           = &typeid ( int );
  types["uint"]          = &typeid ( unsigned int );
  types["double"]        = &typeid ( double );
  types["string"]        = &typeid ( string );
  types["bool"]          = &typeid ( bool );

  string line;

  while ( getline ( stream, line ) ) {

    std::string::size_type hash = line.find ( "#" ); // Cut comments
    if ( hash != string::npos ) line = line.substr ( 0, hash );

    stringstream temp ( line );
    string type;
    string identifier;
    string value;

    temp >> type;
    if ( type == "" ) continue;

    // Check if type is a known type
    if ( types.find ( type ) == types.end () ) throw std::invalid_argument ( aol::strprintf ( "Read: type %s not recognized in File %s at line %d.", type.c_str(), __FILE__, __LINE__ ).c_str() );

    temp >> identifier;
    getline ( temp, value );

    if ( _values.find ( identifier ) != _values.end () ) throw std::invalid_argument ( aol::strprintf ( "variable %s is not unique. File %s line %d.", identifier.c_str(), __FILE__, __LINE__ ).c_str() );

    _values [identifier] = Variable ( types [type], identifier, string ( value ), type );
  }
  stream.close();
}


PParser::Variable PParser::getVariable ( string identifier ) const {
  std::map<string, Variable>::const_iterator found = _values.find ( identifier );
  if ( found == _values.end () ) throw std::invalid_argument ( aol::strprintf ( "variable %s not found in parameter file in File %s at line %d.", identifier.c_str() ,__FILE__, __LINE__ ).c_str() );
  return ( * ( found ) ).second;
}

void PParser::setVariable ( const std::string& identifier, Variable& variable ) {
  std::map<string, Variable>::iterator found = _values.find ( identifier );
  if ( found == _values.end () ) throw std::invalid_argument ( aol::strprintf ( "variable not found in parameter file in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  _values[identifier] = variable;
}


bool PParser::hasVariable ( string identifier ) const {
  std::map<string, Variable>::const_iterator found = _values.find ( identifier );
  return ( found != _values.end () );
}

void PParser::dump ( string filename ) const {
  ofstream outputFile ( filename.c_str () );
  if ( !outputFile.is_open () ) throw std::invalid_argument ( aol::strprintf ( "failed to open in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  for ( std::map<string, Variable>::const_iterator iter = _values.begin(); iter != _values.end(); ++iter )  {
    string tmp;
    tmp = iter->second.getValueString();
    // Trim leading whitespaces
    tmp.erase ( tmp.begin (), std::find_if ( tmp.begin (), tmp.end (), std::not1 ( std::ptr_fun<int, int> ( std::isspace ) ) ) );
    outputFile << setw(15) << left << iter->second.getTypeString() << setw(30) << left << iter->first << tmp << endl;
  }
  outputFile.close();
}

void PParser::dumpToFile ( const string& FileName, const string& Directory ) const {
  string outFileName;
  outFileName += Directory; outFileName += FileName;
  dump ( outFileName );
}


void PParser::addCounterToSaveDirectory ( const string counterFileName, const string additionalName ) {
  if ( aol::fileExists ( counterFileName ) == false ) {
    std::ofstream out ( counterFileName.c_str() );
    out << 0 << endl;
    out.close ( );
  }
  std::fstream counterFile;
  counterFile.open ( counterFileName.c_str() );
  if ( counterFile.is_open () ) {
    string temp;
    std::getline ( counterFile, temp );
    int counter = atoi ( temp.c_str () );
    counterFile.seekg ( ios::beg );
    counterFile << ++counter;
    
    std::string newFileName;
    this->get<string> ( "saveDirectory", newFileName );
    newFileName += additionalName;
    newFileName += "-";
    newFileName += std::to_string(counter);
    Variable var( &typeid ( string ), "saveDirectory", newFileName, "string" );
    this->setVariable ( "saveDirectory", var );
    
  }
  else throw std::invalid_argument ( aol::strprintf ( "cannot open counter file for writing in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  counterFile.close ();
}


// Class PParser::Variable
PParser::Variable::Variable () : _type ( 0 ), _identifier (), _value (), _typeName () {}

PParser::Variable::Variable ( const PParser::Variable& variable ) : _type ( variable._type ), _identifier ( variable._identifier ), _value ( variable._value ), _typeName ( variable._typeName ) {}

PParser::Variable& PParser::Variable::operator = ( const PParser::Variable& variable ) {
  // Beware of self-assignment
  if ( this == &variable ) return *this;

  _type = variable._type;
  _identifier = variable._identifier;
  _value = variable._value;
  _typeName = variable._typeName;

  return *this;
}

PParser::Variable::~Variable () {}

PParser::Variable::Variable ( TypePtr type, string identifier, string value, string typeName )
    : _type ( type ), _identifier ( identifier ), _value ( value ), _typeName ( typeName ) {}

const string & PParser::Variable::getValueString () const { return _value;}
const string & PParser::Variable::getTypeString () const { return _typeName;}


}//end namespace
