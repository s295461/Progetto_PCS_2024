#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

namespace FractureNetwork {


bool ImportFracture(const string& filePath, DiscreteFractureNetwork& fracture)
{
    cout << "FR3_data" << endl;
    string fileNameFR3 = "/FR3_data.txt";
    if(!readFracture(filePath, fileNameFR3, fracture))
        return false;


    PlaneIntersection(fracture);

    //cout << endl;
    //cout << "FR10_data" << endl;
    //string fileNameFR10 = "/FR10_data.txt";
    //readFracture(filePath, fileNameFR10, fracture);
    return true;
}


// Questa funzione apre un file, legge il contenuto e lo salva in strutture dati adeguate
bool readFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture)
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
//    fracture.vertices.reserve(fracture.numFracture);

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

        // Salvo in una matrice le coordinate dei vertici della frattura
        MatrixXd verticesFracture(3,fracture.NumVertices[i]);
        for(unsigned int n = 0; n < fracture.NumVertices[i]; n++)
        {
            verticesFracture(0, n) = x[n];
            verticesFracture(1, n) = y[n];
            verticesFracture(2, n) = z[n];
        }

        // Salvo la matrice di coordinate in un vettore di matrici
        fracture.vertices.push_back(verticesFracture);
    }

    // Chiudo il file__
    file.close();
    return true;
}


// Funzione che calcola il prodotto vettoriale
Vec3d crossProduct(const Vec3d& a, const Vec3d& b)
{
    Vec3d result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}


// Funzione che calcola la norma 2 di un vettore
double norm2(const Vec3d& a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}


// Questa funzione calcola l'intersezione tra due piani creati da due fratture
bool PlaneIntersection(const DiscreteFractureNetwork fracture)
{

    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        unsigned int n = 0;
        for(unsigned int j = 0; j < fracture.numFracture; j++)
        {
            n = 0;
            if( i == j)
            {
                n++;
            }

            if(n == 0)
            {
                unsigned int Id1 = fracture.fractureID[i];
                unsigned int Id2 = fracture.fractureID[j];
                cout << "ID2: " << Id2 << endl;
                Vec3d P;
                P.x = fracture.vertices[Id1](0,0);
                P.y = fracture.vertices[Id1](1,0);
                P.z = fracture.vertices[Id1](2,0);
                Vec3d Q;
                Q.x = fracture.vertices[Id1](0,1);
                Q.y = fracture.vertices[Id1](1,1);
                Q.z = fracture.vertices[Id1](2,1);
                Vec3d R;
                R.x = fracture.vertices[Id1](0,2);
                R.y = fracture.vertices[Id1](1,2);
                R.z = fracture.vertices[Id1](2,2);
                Vec3d A;
                A.x = fracture.vertices[Id2](0,0);
                A.y = fracture.vertices[Id2](1,0);
                A.z = fracture.vertices[Id2](2,0);
                Vec3d B;
                B.x = fracture.vertices[Id2](0,1);
                B.y = fracture.vertices[Id2](1,1);
                B.z = fracture.vertices[Id2](2,1);
                Vec3d C;
                C.x = fracture.vertices[Id2](0,2);
                C.y = fracture.vertices[Id2](1,2);
                C.z = fracture.vertices[Id2](2,2);
                cout << "P: " << P.x << ", R: " << R.x << endl;
                cout << "Q: " << Q.x << ", R: " << R.x << endl;
                Vec3d u1;
                u1.x = P.x - R.x;
                u1.y = P.y - R.y;
                u1.z = P.z - R.z;
                cout << "u1: " << u1.x << endl;

                Vec3d v1;
                v1.x = Q.x - R.x;
                v1.y = Q.y - R.y;
                v1.z = Q.z - R.z;
                cout << "v1: " << v1.x << endl;

                cout << "A: " << A.x << ", C: " << C.x << endl;
                cout << "B: " << B.x << ", C: " << C.x << endl;
                Vec3d u2;
                u2.x = A.x - C.x;
                u2.y = A.y - C.y;
                u2.z = A.z - C.z;
                cout << "u2: " << u2.x << endl;

                Vec3d v2;
                v2.x = B.x - C.x;
                v2.y = B.y - C.y;
                v2.z = B.z - C.z;
                cout << "v2: " << v2.x << endl;

                Vec3d n1;
                n1.x = crossProduct(u1, v1).x / (norm2(u1) * norm2(v1));
                n1.y = crossProduct(u1, v1).y / (norm2(u1) * norm2(v1));
                n1.z = crossProduct(u1, v1).z / (norm2(u1) * norm2(v1));
                cout << "n1.x: " << n1.x << endl;
                cout << "n1.y: " << n1.y << endl;
                cout << "n1.z: " << n1.z << endl;

                Vec3d n2;
                n2.x = crossProduct(u2, v2).x / (norm2(u2) * norm2(v2));
                n2.y = crossProduct(u2, v2).y / (norm2(u2) * norm2(v2));
                n2.z = crossProduct(u2, v2).z / (norm2(u2) * norm2(v2));
                cout << "n2.x: " << n2.x << endl;
                cout << "n2.y: " << n2.y << endl;
                cout << "n2.z: " << n2.z << endl;

                Vec3d t;
                t.x = crossProduct(n1, n2).x;
                t.y = crossProduct(n1, n2).y;
                t.z = crossProduct(n1, n2).z;

                double d1;
                d1 = n1.x * R.x + n1.y * R.y + n2.z * R.z;
                cout << "d1: " << d1 << endl;

                double d2;
                d2 = n2.x * C.x + n2.y * C.y + n2.z * C.z;
                cout << "d2: " << d2 << endl;

                MatrixXd coeff(3, 3);
                coeff(0,0) = n1.x;
                coeff(0,1) = n1.y;
                coeff(0,2) = n1.z;
                coeff(1,0) = n2.x;
                coeff(1,1) = n2.y;
                coeff(1,2) = n2.z;
                coeff(2,0) = t.x;
                coeff(2,1) = t.y;
                coeff(2,2) = t.z;

                Vector3d term;
                term(0) = d1;
                term(1) = d2;
                term(2) = 0;

                Vector3d solution;
                solution = coeff.partialPivLu().solve(term);

            }



        }
    }
    return true;
}

//bool esisteIntersezione(fracture1, fracture2)

// if(!esisteIntersezione)

}
