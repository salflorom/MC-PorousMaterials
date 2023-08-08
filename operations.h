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

			for (int i=1; i<pow; i++) result *= value;
			return result;
		}
		// Function extracted from url "https://cplusplus.com/forum/general/255896"
		// Consultation date: 08/07/2023.
		// Convergence restrictions: abs(x)<1 and c not a negative integer or zero.
		double Hyp2F1(double a, double b, double c, double x){
			const double TOLERANCE = 1e-10;
			double term = a * b * x / c;
			double value = 1.0 + term;
			int n = 1;

			while ( abs( term ) > TOLERANCE ){
				a++, b++, c++, n++;
				term *= a * b * x / c / n;
				value += term;
			}
			return value;
		}
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _OPERATIONS_H_ */

