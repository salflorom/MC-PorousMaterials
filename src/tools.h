#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <string> // string, getline, c_str
#include <cstring> // strlen
#include <sstream> // stringstream
#include "Eigen/Dense"

#include "MC.h"

using namespace std;

/* **************************************************************************** */
class Tools {
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		Tools(void){}
		string LowerCase(string line){
			char lowerLine[50];
			size_t length;
			string str(line);

			length = strlen(line.c_str());
			length = str.copy(lowerLine, length);
			lowerLine[length] = '\0';
			for (size_t i=0; i<length; i++){lowerLine[i] = tolower(line[i]);}
			if (lowerLine[length-1] == ';') lowerLine[length-1] = '\0';
			return lowerLine;
		}
		Eigen::Matrix<string,50,1> SplitString(string str, char delimiter){
			string word;
			Eigen::Matrix<string,50,1> strArray;
			int i=0;
			stringstream ss(str);

			while (!ss.eof()){
				getline(ss, word, delimiter);
				strArray(i) = word;
				i++;
			}
			return strArray;
		}
		template <typename type>
		int FindIndex(const Eigen::MatrixBase<type>& array, int length, string seek){
			for (int i=0; i<length; i++){
				if (array(i).name == seek) {
					return i;
				}
			}
			return -1;
		}
		double Pow(double value, int pow){
			double result=1.;

			for (int i=1; i<=pow; i++) result *= value;
			return result;
		}
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _TOOLS_H_ */

