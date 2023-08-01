#ifndef _OPERATIONS_H_
#define _OPERATIONS_H_

#include <cstring>
#include <fstream>

using namespace std;

/* **************************************************************************** */
class Operations {
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		Operations(void){}
		string LowerCase(string line){
			char lowerLine[50];
			int length;
			string str(line);

			length = strlen(line.c_str());
			length = str.copy(lowerLine, length);
			lowerLine[length] = '\0';
			for (int i=0; i<length; i++){lowerLine[i] = tolower(line[i]);}
			return lowerLine;
		}
		bool FileExists(string fileName){
			bool exists;

			ifstream file(fileName);
			exists = (file ? true : false);
			file.close();
			return exists;
		}
		double Pow(double value, int pow){
			double result=value;

			for (int i=0; i<pow; i++) result *= value;
			return result;
		}
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _OPERATIONS_H_ */

