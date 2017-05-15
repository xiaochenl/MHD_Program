/*
 * Matrix.h
 *
 *  Created on: Mar 11, 2016
 *      Author: majun
 */
#pragma once

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdlib.h>
#include <cmath>
#include "string.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

template<typename Type>

class matrix
{
public:
	matrix(const int row = 1,const int column = 1);
	matrix(const matrix<Type>& rhs);

	~matrix(void);

	void Resize(int row,int column);

	inline Type* & operator [] (int row) const
	{
		return m_RowStart[row];
	}

	void Show(void);

	//Write matrix into file 'Name + N'
	void WriteToFile(const string& fileNameFixed, const int fileNameIndex);

	//Write matrix into file 'Name'
	void WriteToFile(const string& fileNameFixed);
public:
	Type * m_FirstAdd;
	Type ** m_RowStart;

	int m_Row,m_Column;
};


template<typename Type>
matrix<Type>::matrix(const int row,const int column)
{
	m_Row = row;
	m_Column = column;

	m_FirstAdd = new Type [row*column];
	m_RowStart = new Type*[row];

	if (NULL==m_FirstAdd||NULL==m_RowStart) {
		cout << "???ERROR: Memory Exausted!!!" << endl;
		exit(1);
	}

	for (int i=0; i<row; ++i) {
		m_RowStart[i] = m_FirstAdd + i*column;
	}
}


template<typename Type>
matrix<Type>::matrix(const matrix<Type>& rhs)
{
	m_Row = rhs.m_Row;
	m_Column = rhs.m_Column;

	m_FirstAdd = new Type [m_Row*m_Column];
	m_RowStart = new Type*[m_Row];

	if (NULL==m_FirstAdd||NULL==m_RowStart) {
		cout << "???ERROR: Memory Exausted!!!" << endl;
		exit(1);
	}

	for (int i=0; i<m_Row; ++i) {
		m_RowStart[i] = m_FirstAdd + i*m_Column;
	}
}


template<typename Type>
matrix<Type>::~matrix(void)
{
	delete [] m_RowStart;
	m_RowStart = NULL;

	delete [] m_FirstAdd;
	m_FirstAdd = NULL;
}


template<typename Type>
void matrix<Type>::Resize(int row,int column)
{
	delete [] m_RowStart;
	delete [] m_FirstAdd;

	m_Row = row;
	m_Column = column;

	m_FirstAdd = new Type [m_Row*m_Column];
	m_RowStart = new Type*[m_Row];

	if (NULL==m_FirstAdd||NULL==m_RowStart) {
		cout << "???ERROR: Memory Exausted!!!" << endl;
		exit(1);
	}

	for (int i=0; i<m_Row; ++i) {
		m_RowStart[i] = m_FirstAdd + i*m_Column;
	}
}



template<typename Type>
void matrix<Type>::Show(void)
{
	for (int row=0; row<m_Row; ++row) {
		for (int col = 0; col<m_Column; ++col) {
			cout << m_RowStart[row][col] << " ";
		}
		cout << endl;
	}
}


template<typename Type>
void matrix<Type>::WriteToFile(const string& fileNameFixed, const int fileNameIndex)
{
	ofstream myFile;

	char c[512];
	string a2 = fileNameFixed;//,b2=".txt";
	//myitoa(fileNameIndex,c,10);
	sprintf(c,"%d",fileNameIndex);
	a2 += c;
	//a2 += b2;

	//myFile.open(fileNameFixed);
	myFile.open(a2.c_str());


	for (int _I = 0; _I < m_Row; _I++)
	{
		for (int _J = 0; _J < m_Column; _J++)
		{
			myFile <<setprecision(16)<< m_RowStart[_I][_J] << " ";
		}
		myFile << endl;
	}

	myFile << endl;
	myFile.close();
}




template<typename Type>
void matrix<Type>::WriteToFile(const string& fileNameFixed)
{
	ofstream myFile;

	myFile.open(fileNameFixed.c_str());


	for (int _I = 0; _I < m_Row; _I++)
	{
		for (int _J = 0; _J < m_Column; _J++)
		{
			myFile <<setprecision(16)<< m_RowStart[_I][_J] << " ";
			//myFile <<setprecision(16)<< m_FirstAddr[_I][_J] << " ";
		}
		myFile << endl;
	}

	myFile << endl;
	myFile.close();
}
#endif /* MATRIX_H_ */
