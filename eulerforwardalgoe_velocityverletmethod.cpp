
/*
     Earth-Sun system, Sun (0,0), Earth (1,0)  distance unit: AU
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace arma;
// output file as global variable
ofstream ofile;
// function declarations
void output( double, double, double, double, double);

int main(int argc, char* argv[])
{
//  declarations of variables
    char *outfilename;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 2 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file on same line" << endl;
    //    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  string s(outfilename);
  s+=".txt";
  outfilename = const_cast<char*>(s.c_str());

  ofile.open(outfilename);

  int NumberofSteps = 20000;
  double yr=365.0;
  double FinalTime = 10;
  double Step = FinalTime/((double) NumberofSteps);
  double time = 0.0;

  // planet information

  // Initial values  x = 1.0 AU and vy = 2*pi
  double pi = acos(-1.0);
  double FourPi2 = 4*pi*pi;
  double x =  1.0,y =  0.0;
  double vx = (-5.487608572692462E-03)*yr; double vy = (1.633769725891812E-02)*yr;
  double r = sqrt(x*x+y*y);
  double ax=-FourPi2*x/pow(r,3);
  double ay=-FourPi2*y/pow(r,3);

  //check flag
  //-t: check time usage
  //-e: Euler's algo
  //-v: Velocity Verlet method

  for (int i = 1; i < argc; i++) {
      if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("t") != string::npos)) {
           //Time flag (-t):
//          time( n, h,x,u, v,f,v_LU);
//          solved = 1;
      }
      if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("e") != string::npos)) {
          // -e
          // solving the differential equations using Euler's method
          // x(t+dt)=x(t)+dt*v(t);v(t+dt)=v(t)+dt*a(t)

          // v(t+dt)=v(t)-dt*4*pi*pi*x(t)/r(t)^3 (force direction along the opposite of x-axis);

            while (time <= FinalTime){
              x += Step*vx;
              y += Step*vy;
              vx += Step*ax;
              vy += Step*ay;
              r = sqrt(x*x+y*y);
              ax=-FourPi2*x/pow(r,3);
              ay=-FourPi2*y/pow(r,3);

              time += Step;
              output(time, x, y, vx, vy);   // write to file
            }
            ofile.close();  // close output file

          }
      if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("v") != string::npos)){
          //-v
          // velocity verlet method
          // x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
          // v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];
          double ax_new=0, ay_new=0;
          while (time <= FinalTime){
            x += Step*vx+0.5*Step*Step*ax;
            y += Step*vy+0.5*Step*Step*ay;
            r = sqrt(x*x+y*y);
            ax_new=-FourPi2*x/pow(r,3);
            ay_new=-FourPi2*y/pow(r,3);
            vx += 0.5*Step*(ax+ax_new);
            vy += 0.5*Step*(ay+ay_new);
            ax=ax_new;
            ay=ay_new;
            time += Step;
            output(time, x, y, vx, vy);   // write to file
          }
          ofile.close();  // close output file

          }
    }
           return 0;
  }



//    function to write out the final results
void output(double time, double x, double y, double vx, double vy)
{
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << time;
  ofile << setw(15) << setprecision(8) << x;
  ofile << setw(15) << setprecision(8) << y;
  ofile << setw(15) << setprecision(8) << vx;
  ofile << setw(15) << setprecision(8) << vy << endl;
}  // end of function output
