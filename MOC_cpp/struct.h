/* This file contains the structures necessary for the MOC
   code*/

#include <iostream>
#include <cmath>
#include <string>
#include <assert.h>

using namespace std;

class Coord{
	private:
		int *x;
		int *y;

	public:
		// THE BIG THREE

		// Default Ctor:
		Coord() : x(new int(0)), y(new int(0)) {
			cout << "Coord Default Ctor runs" << endl;
		}

		// Initializing Ctor:
		Coord(int x_in, int y_in) : x(new int(x_in)), y(new int(y_in)) {
			cout << "Coord Initializing Ctor runs" << endl;
		}

		// Copy Ctor:
		Coord(const Coord &coord) {
			x = new int(0);
			y = new int(0);

			*x = *(coord.x);
			*y = *(coord.y);

			cout << "Coord Copy Ctor runs" << endl;
		}

		// Dtor:
		~Coord() {
			delete x;
			delete y;
			x = 0;
			y = 0;

			cout << "Coord Dtor runs" << endl;
		}

		// OPERATORS

		// Assignment Operator:
		Coord &operator= (const Coord &rhs) {
			*x = *(rhs.x);
			*y = *(rhs.y);

			cout << "Coord Assignment Operator Ctor runs" << endl;

			return *this;
		}
};

class Quad {
	private:
		Coord coord[2][2];
		int *height;
		int *length;

	public:
		// THE BIG THREE

		// Default Ctor:
		// note: starts off with 4 subcells (2 X 2)
		Quad() : height(new int(1)), length(new int(1)) {
			for (int i = 1; i < 3; ++i) {
				for (int j = 1; j < 3; ++j) {
					coord[i - 1][j - 1] = Coord(i,j);
				}
			}

			cout << "Quad Default Ctor runs" << endl;
		}

		// Initializing Ctor:
		// note: all cells have the same dimension
		Quad(int h_in, int l_in) : height(new int(h_in)), length(new int(l_in)) {
			for (int i = 1; i < 3; ++i) {
				for (int j = 1; j < 3; ++j) {
					coord[i - 1][j - 1] = Coord(i,j);
				}
			}

			cout << "Quad Initializing Ctor runs" << endl;
		}

		// Dtor:
		~Quad() {
			delete height;
			delete length;
			height = 0;
			length = 0;

			cout << "Quad Dtor runs" << endl;
		}

};

class Hex {
	private:
		Coord coord[4][6];
		int *height;
		int *length;

	public:
		// THE BIG THREE

		// Default Ctor:
		// note: starts off with 24 subcells (4 X 6)
		Hex() : height(new int(1)), length(new int(1)) {
			for (int i = 1; i < 5; ++i) {
				for (int j = 1; j < 7; ++j) {
					coord[i - 1][j - 1] = Coord(i,j);
				}			
			}

			cout << "Hex Default Ctor runs" << endl;
		}

		// Initializing Ctor:
		Hex(int h_in, int l_in) : height(new int(h_in)), length(new int(l_in)) {
			for (int i = 1; i < 5; ++i) {
				for (int j = 1; j < 7; ++j) {
					coord[i - 1][j - 1] = Coord(i,j);
				}			
			}

			cout << "Hex Initializing Ctor runs" << endl;
		}

		// Dtor:
		~Hex() {
			delete height;
			delete length;
			height = 0;
			length = 0;

			cout << "Hex Dtor runs" << endl;
		}
};

