#include<math.h>
//Romans class to produce random numbers
class RandomNumberGenerator {
private:
 float width;
 float x1, x2, w;
 float y[2];
public:
 RandomNumberGenerator() { y[0]=0.0; y[1]=0.0;}
 float* Gauss(const float width) {
   //random number generation a la
   //Dr. Everett (Skip) F. Carter Jr.
   //http://www.taygeta.com/random/gaussian.html
      for (int irand = 0; irand < 100; irand++){
     do {
       // generate pair of random numbers
       // generate first random number
       float frand = rand();
       frand = frand/RAND_MAX;
       x1 = 2.0*(frand) - 1.0;
       // generate second random number
       frand = rand();
       frand = frand/RAND_MAX;
       x2 = 2.0*frand - 1.0;
       w = x1 *x1 + x2*x2;
     } while ( w >= 1.0 );
     w = sqrt( (-2.0 * log( w ) ) / w );
     y[0] = width*x1 * w + 1.0;
     y[1] = width*x2 * w + 1.0;
     //cout << width << endl;
     // random->Fill(y1);
     //    random->Fill(y2);
   }
   return &y[0];
 }
 float* EqualDistribution(const float width){
   y[0] = -1000.;
      //assign dummy value to y[1]
   y[1] = -1.;
   while ( y[0] < (1.-width) || y[0] > (1.+width) )
     {
         y[0] = 0.5+rand()/(RAND_MAX*1.0);
     }
   //cout << y[0] << endl;
   //y[0] = width* rand()/RAND_MAX;
   return &y[0];
 }

};
