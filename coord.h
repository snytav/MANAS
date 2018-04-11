/*
 * coord.h
 *
 *  Created on: May 28, 2016
 *      Author: snytav
 */

#ifndef _H_
#define COORD_H_

template<typename T,int extent,int unity>
class Coord{
  T Z;
  public:

     T x,y;
     T z(){return ((extent == 3)*Z + (extent == 2)*((T)unity));}
     void setZ(T zz){Z = zz;}
};

#endif /* COORD_H_ */
