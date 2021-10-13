#ifndef VOBJECT_
#define VOBJECT_ 1

#include <string>

class VObject{
private:
  int vn{};                     // variable number
  std::string *varNames{};      // array of variable name
  double *varValues{};          // array of values of variables

public:
  VObject(const VObject&) = delete;
  VObject& operator=(const VObject&) = delete;
  VObject(int NoOfVars) ;      // constructor : creates arrays, takes number of variables
  ~VObject() ;                 // destructor  : deletes arrays

  // set methods return 0 on succes and -1 on error (+ error message)
  int SetName(int VarNo, std::string VarName) ;    // set variable name no. 'VarNo' to 'VarName'
  int SetNames(std::string *Names) ;               // as above at once (not recommended to use, no error handling)
  
  int SetValue(int VarNo, double Value) ;         // set variable value no. 'VarNo' to 'Value'
  int SetValue(std::string VarName, double Value); // overloaded: set variable 'VarName' to value 'Value'

  // get methods return the value on succes and 0 on error (+ error message!)
  double GetValue(int VarNo) const;               // returns value of variable no. 'VarNo'
  double GetValue(std::string VarName) const;      // overloaded: returns value of variable 'VarName'

  int GetNoOfVariables() const;                   // returns vn (number of variables)
  std::string GetName(int VarNo) const;            // returns the name of variable no. 'VarNo' of "" + error
  std::string* GetNames() const;                   // returns all variable names in an array of size vn

};

#endif
