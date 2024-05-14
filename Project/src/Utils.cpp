#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

namespace FractureNetwork {

bool ImportFracture(const string& filePath, DiscreteFractureNetwork fracture)
{
    cout << "FR3_data" << endl;
    string fileNameFR3 = "/FR3_data.txt";
    readFracture(filePath, fileNameFR3, fracture);

    //cout << endl;
    //cout << "FR10_data" << endl;
    //string fileNameFR10 = "/FR10_data.txt";
    //readFracture(filePath, fileNameFR10, fracture);

    return true;
}


// Questa funzione apre un file, legge il contenuto e lo salva in strutture dati adeguate
bool readFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork fracture)
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
        size_t found = line.find(limit); // found rappresenta la posizione del ";"
        while(found != string::npos) //Cerco finchè non arrivo alla fine della riga
        {
            line.erase(found, 1); // Cancello il ";"
            found = line.find(limit, found); // Cerco il prossimo ";"
        }
        listLines.push_back(line);
    }
    // Elimino la prima riga
    listLines.pop_front();

    // Verifico che il file non sia vuoto
    unsigned int dimension = listLines.size();
    if(dimension == 0)
    {
        cerr << "The discrete fracture network does not exist" << endl;
        return false;
    }

    // Salvo il numero di fratture
    istringstream converterNum(listLines.front());
    converterNum >> fracture.numFracture;

    // Elimino la riga contenente il numero di fratture perchè non mi serve più
    listLines.pop_front();

    // Inizializzo la memoria
    fracture.fractureID.reserve(fracture.numFracture);
    fracture.NumVertices.reserve(fracture.numFracture);
    fracture.vertices.reserve(fracture.numFracture);

    // Ciclo sul numero di fratture nel file
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        // Elimino la prima riga che non mi serve
        listLines.pop_front();

        // Memorizzo id e numero di vertici
        unsigned int id;
        unsigned int num;
        istringstream converterID(listLines.front());
        converterID >> id >> num;
        cout << "id: " << id << endl;
        cout << "num: " << num << endl;


        fracture.fractureID.push_back(id);
        fracture.NumVertices.push_back(num);

        cout << "ID[" << i << "]:" << fracture.fractureID[i] << endl;
        cout << "NUM[" << i << "]:" << fracture.NumVertices[i] << endl;

        // Elimino la riga di id e numero di vertici e la successiva che non mi serve
        listLines.pop_front();
        listLines.pop_front();

        VectorXd x(fracture.NumVertices[i]);
        VectorXd y(fracture.NumVertices[i]);
        VectorXd z(fracture.NumVertices[i]);

        // Salvo le coordinate x dei vertici
        istringstream converterCoordX(listLines.front());
        for(unsigned int n = 0; n < fracture.NumVertices[i]; n++)
        {
            converterCoordX >> x[n];

        }
        listLines.pop_front(); // Elimino la riga delle x

        // Salvo le coordinate y dei vertici
        istringstream converterCoordY(listLines.front());
        for(unsigned int n = 0; n < fracture.NumVertices[i]; n++)
        {
            converterCoordY >> y[n];

        }
        listLines.pop_front(); // Elimino la riga delle y

        // Salvo le coordinate z dei vertici
        istringstream converterCoordZ(listLines.front());
        for(unsigned int n = 0; n < fracture.NumVertices[i]; n++)
        {
            converterCoordZ >> z[n];

        }
        listLines.pop_front(); // Elimino la riga delle z

        // Stampo i vertici
        for(unsigned int n = 0; n < fracture.NumVertices[i]; n++)
        {
            cout << scientific << setprecision(16) << "Vertices " << n+1 << ": x: " << x[n] << "; y: " << y[n] << "; z: " << z[n] << endl;
        }

        // Salvo in memoria i verici
        vector<VectorXd> coord;
        coord.push_back(x);
        coord.push_back(y);
        coord.push_back(z);
        fracture.vertices.push_back(coord);

        // Salvarli in una matrice??

        // cout << fracture.vertices << endl;

        //cout << scientific << setprecision(16) << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
        //cout << scientific << setprecision(16) << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;
        //cout << scientific << setprecision(16) << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << endl;

    }

    file.close();
    return true;

}




}
