#include <map>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include "PMT.h"
#include "PMTParam.h"
// DD4hep include files
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "Math/Vector3D.h"
#include "Math/Transform3D.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <regex>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;


void trim(string &s)
{
  if( !s.empty() )
  {
    s.erase(0,s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
  }
}

// split string by delim, one example is
// a=s_split("adfasd,fasd,fa ,,asd,fasdf,asdf   , ,", "[ ,]+")
std::vector<std::string> s_split( std::string &in, const std::string &delim)
{
  std::regex re{delim};
  trim(in);
  return std::vector<std::string>{
      std::sregex_token_iterator(in.begin(), in.end(), re, -1),
      std::sregex_token_iterator()};
}

PMTParam::PMTParam(std::string positionFilePath, std::string typeFilePath)
{
  this->_positionFilePath = positionFilePath;
  this->_typeFilePath = typeFilePath;
  readFileGetPosition(6);
  readFileGetType(2);

}

PMTParam::PMTParam(std::string positionFilePath, PMTType defaultPMTType)
{
  this->_positionFilePath = positionFilePath;
  this->_defaultPMTType = defaultPMTType;
  readFileGetPosition(6);
}

PMTParam::~PMTParam()
{
}

/// @brief
/// @param paramNumPerLine the num of params split by space or tab per line in the file
/// @return void. but save map<int, PMT> to  _pmtParam
void PMTParam::readFileGetPosition(int paramNumPerLine)
{
  // vector<PMT> pmtPos;
  string filePath=_positionFilePath;
  std::ifstream input(filePath.c_str());
  std::string tmp_line;
  while (input.good())
  {
    std::getline(input, tmp_line);
    if (tmp_line == "") continue;

    vector<string> num = s_split(tmp_line, "[ \t]+");
    if(num.size() != uint64_t(paramNumPerLine))
    {
      cout << "@@@ PMTParam " << __LINE__ << ", file=" << filePath <<  ", temp_line=\'" << tmp_line << "\', num size=" << num.size() << endl;
      continue;
    }
    else
    {
      int idx = stoi(num[0]);
      Position pos(stof(num[1]), stof(num[2]), stof(num[3]));
      XYZAngles ang(stof(num[4]), stof(num[5]), 0);
      _pmtParam[idx] = PMT(pos, ang, _defaultPMTType);
    }
  }
  cout << "@@@ PMTParam " << __LINE__ <<  ", file=" << filePath << ", line nums=" << _pmtParam.size() << endl;
  // return pmtPos;
}

void PMTParam::readFileGetType(int paramNumPerLine)
{
  string filePath=_typeFilePath;
  std::ifstream input(filePath.c_str());
  std::string tmp_line;

  while (input.good())
  {
    std::getline(input, tmp_line);
    if (tmp_line == "") continue;

    vector<string> num = s_split(tmp_line, "[ \t]+");
    if(num.size() != uint64_t(paramNumPerLine))
    {
      cout << "@@@ PMTParam " << __LINE__ << ", file=" << filePath << ", temp_line=\'" << tmp_line << "\', num size=" << num.size() << endl;
      continue;
    }
    else
    {
      int idx = stoi(num[0]);
      PMTType tmp;
      if(num[1] == "Hamamatsu")
        tmp = PMTType::Hamamatsu;
      else if(num[1] == "HighQENNVT")
        tmp = PMTType::HighQENNVT;
      else if(num[1] == "NNVT")
        tmp = PMTType::NNVT;
      else
        tmp = PMTType::NotSet;
      _pmtParam[idx].setPMTType(tmp);
    }
  }
  cout << "@@@ PMTParam " << __LINE__ << ", file=" << filePath << ", line nums=" << _pmtParam.size() << endl;
  // return pmtPos;
}

map<int,PMT> PMTParam::getPMTParam()
{
  return _pmtParam;
}
