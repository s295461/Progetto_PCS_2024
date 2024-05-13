#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string& filePath)
{
    string fileNameFR3 = "/FR3_data.txt";
    readFracture(filePath, fileNameFR3);


    string fileNameFR10 = "/FR10_data.txt";
    readFracture(filePath, fileNameFR10);

    return true;
}


bool readFracture(const string& filePath, const string& fileName)
{
    // Apro il file
    ifstream file;
    file.open(filePath + fileName);

    // Verifico che si apra correttamente
    if(file.fail())
    {
        cout << "Err" << endl;
        return false;

    }

    // Leggo il file e sostituisce "; " con " ", inoltre salvo il contenuto in una lista di stringhe
    list<string> listLines;
    string line;
    string limit = "; ";
    while(getline(file, line))
    {
        cout << line << endl;
        size_t found = line.find(limit); // found rappresenta la posizione del ";"
        while(found != string::npos) //Cerco finchè non arrivo alla fine della riga
        {
            line.replace(found, 1, " "); // Sostituisco ";" con " "
            found = line.find(limit, found); // Cerco il prossimo ";"
        }
        cout << line << endl;
        listLines.push_back(line);
    }
    // Elimino la prima riga
    listLines.pop_front();

    file.close();
    return true;

}




}
