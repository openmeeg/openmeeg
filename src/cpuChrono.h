#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

class cpuChrono
{
private:
    clock_t ellapsed;
    clock_t tstart;
public:
    cpuChrono(){ellapsed=0; tstart=0;}
    ~cpuChrono(){}
    void start ()
    {
        tstart = clock();
    }
    void stop()
    {
        ellapsed += (ellapsed+clock()-tstart);
    }
    void zero()
    {
        ellapsed=0;
    }
    clock_t getEllapsedT ()
    {
        return ellapsed;
    }
    double getEllapsedS ()
    {
        return (double)(ellapsed) / CLOCKS_PER_SEC;
    }
    void dispEllapsed ()
    {
        cout <<  "------------------------------------------" << endl;
        cout <<  "| Elapsed Time: " << getEllapsedS() << " s." << endl;
        cout <<  "------------------------------------------" << endl;
    }
};
