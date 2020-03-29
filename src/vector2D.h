#ifndef __vector2D_h__
#define __vector2D_h__
 
 class vector2D
 {
	public:
		
		vector2D(void);										// Default constructor
		vector2D(float,float);								// Constructor via component
			
		float &operator[](int);								// For array-style access
		vector2D operator*(float scale);					// For scalar multiplication
		
		float x;											
		float y;
};

#endif

