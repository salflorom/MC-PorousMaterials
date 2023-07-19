#ifndef _OPERATIONS_H_
#define _OPERATIONS_H_

#include <cstring>

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
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _OPERATIONS_H_ */

