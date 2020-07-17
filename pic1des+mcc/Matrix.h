#ifndef Matrix_h
#define Matrix_h

//		Class		: Matrix (template version)
//		Version		: 2.4.1
//		Programmer	: Jun Hasegawa 
//		Last update	: December 23, 2006
//
//		Copyright (c) 2006-2013 Tokyo Tech. All rights reserved.


#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory>
using namespace std;

const long max_matrix_size = 10000; // upper limit

template<class T>
class Matrix {
public:
	Matrix() throw() : nr_(0), ptr_(NULL) {}	// default constructor
	Matrix( long nr, long nc = 0 ) throw( bad_alloc, length_error );
	Matrix( long nr, long nc, const T& a ) throw( bad_alloc, length_error );	//_
	Matrix( const Matrix& m ) throw( bad_alloc );
	virtual ~Matrix() throw(); 

	// operator[] for fast data access without range check
	vector<T>& operator[]( long i ) throw() { return ptr_[i]; }
	const vector<T>& operator[]( long i ) const throw() { return ptr_[i]; }

	// operator =, +=, -=
	Matrix& operator=( const Matrix& m ) throw();
	Matrix<T>& operator+=( const Matrix& m ) throw( bad_alloc, logic_error );
	Matrix<T>& operator-=( const Matrix& m ) throw( bad_alloc, logic_error );
	Matrix<T>& operator+=( const T& a ) throw( bad_alloc );
	Matrix<T>& operator-=( const T& a ) throw( bad_alloc );

	// member functions for data access with range check
	T& at( long i, long j ) throw( out_of_range );
	const T& at( long i, long j ) const throw( out_of_range );
	vector<T>& at( long i ) throw( out_of_range );
	const vector<T>& at( long i ) const throw( out_of_range );

	long nr() const throw() { return nr_; }	// return number of rows
	long nc( long i = 0 ) const throw() { return at(i).size(); }	// return number of colomns
	void push_back( const vector<T>& a );	// push back data array to the end of the colomns
	void clear() throw();	// delete all data from the matrix
	void fill( const T& a ) throw();	// fill all data with a
	void magnify( const T& a ) throw();	// magnify all data with a
	void resize( long nr, long nc = 0 ) throw( bad_alloc, length_error );	// resize the matrix
	void transpose() throw( bad_alloc );	// transpose the matrix
	bool empty() { return ( ptr_ == NULL ); }	// check if the matrix is empty

	// member functions for data output
	void show() const throw();	// show all data to console
	void output( ostream& out ) const throw(); // output all data to stream
	void input( istream& in, long nr, long nc = 0 ); // input data from stream

	void show_transposed() const throw();	// show all transposed data to console
	void output_transposed( ostream& out ) const throw(); // output all transposed data to stream


protected:
	vector<T>* init( long nr, long nc = 0 ) throw( bad_alloc, length_error );
	void copy( const Matrix& m ) throw( bad_alloc );

	long nr_;	// number of lines
	vector<T>* ptr_;
};

template<class T>
Matrix<T>::Matrix( long nr, long nc ) throw( bad_alloc, length_error )
: nr_(nr), ptr_(init(nr,nc))
{
	//
}

template<class T>
Matrix<T>::Matrix( long nr, long nc, const T& a ) throw( bad_alloc, length_error )
: nr_(nr), ptr_(init(nr,nc))
{
	fill( a );
}

template<class T>
Matrix<T>::Matrix( const Matrix& m ) throw( bad_alloc )
: nr_(0), ptr_(NULL)
{
	copy( m );
}

template<class T>
Matrix<T>::~Matrix() throw()
{
	if ( ptr_ != NULL )
		delete[] ptr_;
}

template<class T>
Matrix<T>& Matrix<T>::operator=( const Matrix& m ) throw()
{
	if ( this == &m ) return *this;	// check self substitution
	copy( m );
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=( const Matrix& m ) throw( bad_alloc, logic_error )
{
	if ( this == &m ) return *this;	// check self substitution
	if ( ptr_ == NULL ) copy( m );
	else {
		if ( nr_ != m.nr() )
			throw logic_error( "type mismatch in Matrix<T>::operator+=()" );
		for ( long i = 0; i < nr_; i++ ) {
			if ( nc(i) != m.nc(i) ) throw logic_error( "type mismatch in Matrix<T>::operator+=()" );
			for ( long j = 0; j < nc(i); j++ )
				ptr_[i][j] += m.ptr_[i][j];
		}
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=( const T& a ) throw( bad_alloc )
{
	if ( ptr_ == NULL ) {
		nr_ = 1;
		ptr_ = init( 1 );
		ptr_[0].push_back( a );
	}
	else {
		for ( long i = 0; i < nr_; i++ )
			for ( long j = 0; j < nc(i); j++ )
				ptr_[i][j] += a;
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=( const Matrix& m ) throw( bad_alloc, logic_error )
{
	if ( this == &m ) return *this;	// check self substitution
	if ( ptr_ == NULL ) {
		copy( m );
		magnify( -1 );
	}
	else {
		if ( nr_ != m.nr() )
			throw logic_error( "type mismatch in Matrix<T>::operator+=()" );
		for ( long i = 0; i < nr_; i++ ) {
			if ( nc(i) != m.nc(i) ) throw logic_error( "type mismatch in Matrix<T>::operator-=()" );
			for ( long j = 0; j < nc(i); j++ )
				ptr_[i][j] -= m.ptr_[i][j];
		}
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=( const T& a ) throw( bad_alloc )
{
	if ( ptr_ == NULL ) {
		nr_ = 1;
		ptr_ = init( 1 );
		ptr_[0].push_back( a );
	}
	else {
		for ( long i = 0; i < nr_; i++ )
			for ( long j = 0; j < nc(i); j++ )
				ptr_[i][j] -= a;
	}
	return *this;
}

template<class T>
T& Matrix<T>::at( long i, long j ) throw( out_of_range )
{
	if ( i < 0 || i > nr_-1 )
		throw out_of_range( "Matrix::at index out of range" ); 
	return ptr_[i].at(j);
}

template<class T>
const T& Matrix<T>::at( long i, long j ) const throw( out_of_range )
{
	if ( i < 0 || i > nr_-1 )
		throw out_of_range( "Matrix::at index out of range" ); 
	return ptr_[i].at(j);
}

template<class T>
vector<T>& Matrix<T>::at( long i ) throw( out_of_range )
{
	if ( i < 0 || i > nr_-1 )
		throw out_of_range( "Matrix::at index out of range" ); 
	return ptr_[i];
}

template<class T>
const vector<T>& Matrix<T>::at( long i ) const throw( out_of_range )
{
	if ( i < 0 || i > nr_-1 )
		throw out_of_range( "Matrix::at index out of range" ); 
	return ptr_[i];
}

template<class T>
void Matrix<T>::push_back( const vector<T>& a )
{
	long n = a.size();
	if ( n > nr_ ) n = nr_;
	for ( long i = 0; i < n; i++ )
		ptr_[i].push_back( a.at(i) );
}

template<class T>
void Matrix<T>::clear() throw()
{
	delete[] ptr_;
	ptr_ = NULL;
	nr_ = 0;
}

template<class T>
void Matrix<T>::fill( const T& a ) throw()
{
	for ( long i = 0; i < nr_; i++ )
		for ( long j = 0; j < ptr_[i].size(); j++ )
			ptr_[i][j] = a;
}

template<class T>
void Matrix<T>::magnify( const T& a ) throw()
{
	for ( long i = 0; i < nr_; i++ )
		for ( long j = 0; j < ptr_[i].size(); j++ )
			ptr_[i][j] *= a;
}

template<class T>
void Matrix<T>::resize( long nr, long nc ) throw( bad_alloc, length_error )
{
	vector<T>* t_ptr = init( nr, nc );
	if ( ptr_ != NULL ) {
		long nr_min = min<long>( nr, nr_ );
		for ( long i = 0; i < nr_min; i++ ) {
			long nc_min = min<long>( nc, ptr_[i].size() );
			if ( nc_min != 0 ) {
				for ( long j = 0; j < nc_min; j++ )
					t_ptr[i][j] = ptr_[i][j];
			}
		}
	}
	delete[] ptr_;
	ptr_ = t_ptr;
	nr_ = nr;
}


template<class T>
void Matrix<T>::transpose() throw( bad_alloc )
{
	if ( ptr_ == NULL ) return;
	long nc_max = 0;
	for ( long i = 0; i < nr_; i++ )
		if ( ptr_[i].size() > nc_max ) nc_max = ptr_[i].size();
	vector<T>* t_ptr = init( nc_max, nr_ );
	for ( long i = 0; i < nr_; i++ )
		for ( long j = 0; j < ptr_[i].size(); j++ )
			t_ptr[j][i] = ptr_[i][j];
	delete[] ptr_;
	ptr_ = t_ptr;	// set new pointer
	nr_ = nc_max;
}


template<class T>
void Matrix<T>::show() const throw()
{
	output( cout );
}

template<class T>
void Matrix<T>::output( ostream& out ) const throw()
{
	if ( ptr_ == NULL ) return;
	for ( long i = 0; i < nr_; i++ ) {
		if ( ptr_[i].size() != 0 ) {
			out << ptr_[i][0];
			for ( long j = 1; j < ptr_[i].size(); j++ ) {
				out << '\t' << ptr_[i][j];
			}
		}
		out << '\n';
	}
}

template<class T>
void Matrix<T>::input( istream& in, long nr, long nc )
{
	nr_ = nr;
	ptr_ = init( nr, nc );
	
	if ( nc == 0 ) {
		for (;;) {
			T buf;
			in >> buf;
			if ( in.eof() ) break;
			ptr_[0].push_back(buf);
			for ( long i = 1; i < nr; i++ ) {
				in >> buf;
				ptr_[i].push_back(buf);
			}
		}
	}
	else {
		for ( long j = 0; j < nc; j++ )
			for ( int i = 0; i < nr; i++ )
				in >> ptr_[i][j];
	}
}

template<class T>
void Matrix<T>::show_transposed() const throw()
{
	output_transposed( cout );
}

template<class T>
void Matrix<T>::output_transposed( ostream& out ) const throw()
{
	if ( ptr_ == NULL ) return;
	long nc_max = 0;
	for ( long i = 0; i < nr_; i++ )
		if ( ptr_[i].size() > nc_max ) nc_max = ptr_[i].size();
	for ( long j = 0; j < nc_max; j++ ) {
		if ( j < ptr_[0].size() ) out << ptr_[0][j];
		for ( long i = 1; i < nr_; i++ ) {
			out << '\t';
			if ( j < ptr_[i].size() ) out << ptr_[i][j];
		}
		out << '\n';
	}
}

// protected member functions

template<class T>
vector<T>* Matrix<T>::init( long nr, long nc ) throw( bad_alloc, length_error )
{
	vector<T>* t_ptr = new vector<T>[nr];
	if ( nc != 0 ) {
		for( long i = 0; i < nr; i++ ) {
			try {
				t_ptr[i].resize( nc );
			}
			catch( bad_alloc ) {
				delete[] t_ptr;
				throw;
			}
			catch( length_error ) {
				delete[] t_ptr;
				throw;
			}
		}
	}
	return t_ptr;
}

template<class T>
void Matrix<T>::copy( const Matrix& m ) throw( bad_alloc )
{
	if ( ptr_ != NULL )
		delete[] ptr_;
	nr_ = m.nr();
	ptr_ = new vector<T>[nr_];
	for ( long i = 0; i < nr_; i++ )
		ptr_[i] = m.ptr_[i];
}

#endif

