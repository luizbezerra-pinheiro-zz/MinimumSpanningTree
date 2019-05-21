// basic file operations
#include <iostream>
#include <fstream>
#include <string>

int main()
{
    std::ofstream weightFile;
    weightFile.open("input/weightFile1.txt");
    srand(time(0));
    for (int i = 0; i < 10000000; i++)
    {
        double w = 1.0*(rand()%100);
        if(w>0)
            weightFile << w << ' ';
    }
    weightFile << '\n';
    weightFile.close();
    return 0;
}