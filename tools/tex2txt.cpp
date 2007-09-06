#include <matrice.h>
#include <symmatrice.h>
#include <vecteur.h>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

int main( int argc, char **argv)
{
    if(argc != 3) 
    {
        cout << "usage: "<< argv[0] <<" input_texture_file output_txt_file" << endl;
        exit(1);
    }

    std::ifstream file(argv[1]);
    if (!file.is_open()) 
    {
        cout << "Problem reading : " << argv[1] << endl;
        return 1;
    }

    vector<string> tokens;

    std::string line;
    while (std::getline(file, line, '\n'))
    {
        if (!line.empty())
        {
            std::istringstream buffer(line);
            string tile;
            while (buffer >> tile)
            {
                tokens.push_back(tile);
            }
        }
    }

    vector<string>::iterator it = tokens.begin();

    if (*it != "ascii")
    {
        cout << "Not ascii format" << endl;
        return 1;
    }
    ++it;

    if (*it != "FLOAT")
    {
        cout << "Not FLOAT format but " << *it << endl;
        return 1;
    }
    ++it;

    int nb_time_steps = atoi((*it).c_str());
    cout << "Nb time steps : " << nb_time_steps << endl;
    ++it ;

    ++it ; // skipping header

    int nb_points = atoi((*it).c_str());
    cout << "Nb points : " << nb_points << endl;

    matrice m(nb_points,nb_time_steps);

    for (int i = 0 ; i < nb_time_steps ; i++)
    {
        nb_points = atoi((*it).c_str());;
        cout << "Nb points at time " << i << " : " << nb_points << endl;
        ++it;
        for (int j = 0 ; j < nb_points ; j++)
        {
            m(j,i) = (float) atof((*it).c_str());
            ++it;
        }
        if (i != (nb_time_steps-1)) ++it;
    }

    m.saveTxt(argv[2]);

    return 0;
}
