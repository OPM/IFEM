#include "Tensor.h"
#include "Vec3.h"
#include "Vec3Oper.h"

/* To compile on Linux:
  g++ -DReal=double -DUSE_CBLAS -I../LinAlg \
  testTensor.cpp Tensor.C Vec3Oper.C -llapack -lblas
*/

int main (int argc, char** argv)
{
  double data[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

  Tensor T2(std::vector<double>(data,data+4)), t2(2);
  Tensor T3(std::vector<double>(data,data+9)), t3(3);

  std::cout <<"T2:\n"<< T2;
  std::cout <<"T3:\n"<< T3;

  t2=T2; std::cout <<"T2.shift(1):\n"<< t2.shift(1);
  t2=T2; std::cout <<"T2.shift(-1):\n"<< t2.shift(-1);
  t2=T2; std::cout <<"T2.shift(2):\n"<< t2.shift(2);
  t2=T2; std::cout <<"T2.shift(-2):\n"<< t2.shift(-2);
  std::cout << std::endl;

  t3=T3; std::cout <<"T3.shift(1):\n"<< t3.shift(1);
  t3=T3; std::cout <<"T3.shift(-1):\n"<< t3.shift(-1);
  t3=T3; std::cout <<"T3.shift(2):\n"<< t3.shift(2);
  t3=T3; std::cout <<"T3.shift(-2):\n"<< t3.shift(-2);
  std::cout << std::endl;

  SymmTensor S2(2);
  SymmTensor S3(3);
  SymmTensor S4(2,true);

  S2(1,1) = 1.0; S2(2,2) = 2.0;
  S3(1,1) = 1.0; S3(2,2) = 2.0; S3(3,3) = 3.0;
  S4(1,1) = 1.0; S4(2,2) = 2.0; S4(3,3) = 0.5;

  std::cout <<"S2:\n"<< S2;
  std::cout <<"S3:\n"<< S3;
  std::cout <<"S4:\n"<< S4;

  const double angle = 25.0*M_PI/180.0;
  T2(1,1) = T2(2,2) = cos(angle);
  T2(1,2) = -sin(angle);
  T2(2,1) =  sin(angle);
  T3 = Tensor(Vec3(1.0,3.0,2.0));
  std::cout <<"T2:\n"<< T2;
  std::cout <<"T3:\n"<< T3;

  S2.transform(T2);
  S3.transform(T3);
  S4.transform(T2);
  std::cout <<"transformed S2:\n"<< S2;
  std::cout <<"transformed S3:\n"<< S3;
  std::cout <<"transformed S4:\n"<< S4;

  Vec3 p;
  S2.principal(p);
  std::cout <<"principal S2: "<< p << std::endl;
  S3.principal(p);
  std::cout <<"principal S3: "<< p << std::endl;
  S4.principal(p);
  std::cout <<"principal S4: "<< p << std::endl;

  Vec3Vec dir;
  S2.principal(p,dir);
  std::cout <<"\nprincipal S2: "<< p;
  std::cout <<"\ndir1: "<< dir[0] <<"\ndir2: "<< dir[1];
  S3.principal(p,dir);
  std::cout <<"\nprincipal S3: "<< p;
  std::cout <<"\ndir1: "<< dir[0] <<"\ndir2: "<< dir[1] <<"\ndir3: "<< dir[2];
  S4.principal(p,dir);
  std::cout <<"\nprincipal S4: "<< p;
  std::cout <<"\ndir1: "<< dir[0] <<"\ndir2: "<< dir[1] <<"\ndir3: "<< dir[2];
  std::cout << std::endl;

  return 0;
}
