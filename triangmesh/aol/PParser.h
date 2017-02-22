#ifndef __PPARSER_H
#define __PPARSER_H

#include <general.h>

namespace aol {

//! A Simple Parameter Parser that can read parameter files of the format "type identifier value".
//! It can easily be extended for new types by adding them to the table in PParser::read
//! Each type must understand the operator >> that is used to interpret the value
//! Run-Time Type safety is guaranteed via typeid ()
//! use as parser.get( "name", variable )
class PParser {

private:
  //! Type information that can be compared
  typedef const type_info* TypePtr;

  /**********************************************************************************/
  //! Storage for a variable type information, identifier, and (still-to-decode) value
  class Variable {

  private:
    //! Type information. Will not be deleted, since memory is owned by run time system
    TypePtr _type;
    //! Name identifier
    string _identifier;
    //! Value as a string, will be decoded to appropriate type by operator >>
    string _value;
    //! Internal type name
    string _typeName;

  public:
    //! Default constructor for STL compatibility
    Variable ();
    //! Copy Constructor
    Variable ( const Variable& variable );
    //! Assignment Operator
    Variable& operator = ( const Variable& variable );
    //! Destructor
    virtual ~Variable ();
    //! Construct from its members
    Variable ( TypePtr type, string identifier, string value, string typeName );

    //! Decode value to appropriate type. Type safety guaranteed at run time by typeid comparison
    template <class DataType> void getValue ( DataType& value ) const {
        if ( *_type != typeid ( DataType ) ) throw std::invalid_argument ( aol::strprintf ( "wrong typeid in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        stringstream temp ( _value );
        temp >> value;
    }
    //! Get value as a string
    const string & getValueString () const;
    //! Get type name as a string
    const string & getTypeString () const;
  };

  //! Get Variable for specified identifier
  Variable getVariable ( string identifier ) const;

public:

  PParser ();
  PParser ( const PParser& parser );
  PParser& operator = ( const PParser& parser );
  virtual ~PParser () ;
  explicit PParser ( const std::string & ParFilename );

  //! Decode value of given identifier to appropriate type. Type safety is guaranteed by Variable::get<DataType>
  template <class DataType> void get ( string identifier, DataType& value ) const {
    Variable variable = getVariable ( identifier );
    variable.getValue ( value );
  }
  template <class DataType> DataType get ( string identifier) const {
    DataType value;
    Variable variable = getVariable ( identifier );
    variable.getValue ( value );
    return value;
  }
  void read ( const std::string& fileName );
  bool hasVariable ( string identifier ) const;
  void setVariable ( const string& identifier, Variable& variable );
  void dump ( string filename ) const;
  void dumpToFile ( const string& FileName, const string& Directory ) const;
  void addCounterToSaveDirectory ( const string, const string additionalName = "" );
private:
  //! Stores Identifier--Variable Pairs
  std::map <string, Variable> _values;
};

}// end namespace

#endif
