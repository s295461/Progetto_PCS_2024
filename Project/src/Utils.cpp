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
// Vec3d crossProduct(const Vec3d& a, const Vec3d& b)
// {
//     Vec3d result;
//     result.x = a.y * b.z - a.z * b.y;
//     result.y = a.z * b.x - a.x * b.z;
//     result.z = a.x * b.y - a.y * b.x;
//     return result;
// }


// Funzione che calcola la norma 2 di un vettore
// double norm2(const Vec3d& a)
// {
//     return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
// }


// Funzione che trova le tracce
//bool findTraces(const Vec3d s, const Vec3d Point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2)
//{

//}

// Questa funzione calcola l'intersezione tra due piani creati da due fratture
bool PlaneIntersection(const DiscreteFractureNetwork fracture)
{

    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        for(unsigned int j = 0; j < fracture.numFracture; j++)
        {
            if(i != j)
            {
                // Seleziono due fratture
                unsigned int Id1 = fracture.fractureID[i];
                unsigned int Id2 = fracture.fractureID[j];

                // Prendo i vertici di tre punti per ogni frattura
                // Vec3d P;
                // P.x = fracture.vertices[Id1](0,0);
                // P.y = fracture.vertices[Id1](1,0);
                // P.z = fracture.vertices[Id1](2,0);
                // Vec3d Q;
                // Q.x = fracture.vertices[Id1](0,1);
                // Q.y = fracture.vertices[Id1](1,1);
                // Q.z = fracture.vertices[Id1](2,1);
                // Vec3d R;
                // R.x = fracture.vertices[Id1](0,2);
                // R.y = fracture.vertices[Id1](1,2);
                // R.z = fracture.vertices[Id1](2,2);
                // Vec3d A;
                // A.x = fracture.vertices[Id2](0,0);
                // A.y = fracture.vertices[Id2](1,0);
                // A.z = fracture.vertices[Id2](2,0);
                // Vec3d B;
                // B.x = fracture.vertices[Id2](0,1);
                // B.y = fracture.vertices[Id2](1,1);
                // B.z = fracture.vertices[Id2](2,1);
                // Vec3d C;
                // C.x = fracture.vertices[Id2](0,2);
                // C.y = fracture.vertices[Id2](1,2);
                // C.z = fracture.vertices[Id2](2,2);

                // Calcolo i vettori
                Vector3d u1;
                u1(0) = fracture.vertices[Id1](0,0) - fracture.vertices[Id1](0,2);
                u1(1) = fracture.vertices[Id1](1,0) - fracture.vertices[Id1](1,2);
                u1(2) = fracture.vertices[Id1](2,0) - fracture.vertices[Id1](2,2);

                Vector3d v1;
                v1(0) = fracture.vertices[Id1](0,1) - fracture.vertices[Id1](0,2);
                v1(1) = fracture.vertices[Id1](1,1) - fracture.vertices[Id1](1,2);
                v1(2) = fracture.vertices[Id1](2,1) - fracture.vertices[Id1](2,2);

                Vector3d u2;
                u2(0) = fracture.vertices[Id2](0,0) - fracture.vertices[Id2](0,2);
                u2(1) = fracture.vertices[Id2](1,0) - fracture.vertices[Id2](1,2);
                u2(2) = fracture.vertices[Id2](2,0) - fracture.vertices[Id2](2,2);

                Vector3d v2;
                v2(0) = fracture.vertices[Id2](0,1) - fracture.vertices[Id2](0,2);
                v2(1) = fracture.vertices[Id2](1,1) - fracture.vertices[Id2](1,2);
                v2(2) = fracture.vertices[Id2](2,1) - fracture.vertices[Id2](2,2);

                // Calcolo le normali ai due piani
                Vector3d n1 = u1.cross(v1) / (u1.norm() * v1.norm());
                //n1(0) = crossProduct(u1, v1).x / (norm2(u1) * norm2(v1));
                //n1(1) = crossProduct(u1, v1).y / (norm2(u1) * norm2(v1));
                //n1(2) = crossProduct(u1, v1).z / (norm2(u1) * norm2(v1));

                Vector3d n2 = u2.cross(v2) / (u2.norm() * v2.norm());
                // n2(0) = u2.cross(v2).x / (norm2(u2) * norm2(v2));
                // n2(1) = crossProduct(u2, v2).y / (norm2(u2) * norm2(v2));
                // n2(2) = crossProduct(u2, v2).z / (norm2(u2) * norm2(v2));

                // Calcolo il vettore direzione della retta di intersezione
                Vector3d s = n1.cross(n2);
                //s.x = crossProduct(n1, n2).x;
                //s.y = crossProduct(n1, n2).y;
                //s.z = crossProduct(n1, n2).z;

                // Calcolo il termine noto
                double d1;
                d1 = n1(0) * fracture.vertices[Id1](0,2) + n1(1) * fracture.vertices[Id1](1,2) + n1(2) * fracture.vertices[Id1](2,2);

                double d2;
                d2 = n2(0) * fracture.vertices[Id2](0,2) + n2(1) * fracture.vertices[Id2](1,2) + n2(2) * fracture.vertices[Id2](2,2);

                // Inserisco i coefficienti del sistema formato dai due piani ricavati dalle fratture e dal piano con normale s in una matrice
                Matrix3d coeff;
                coeff << n1, n2, s;
                coeff(0,0) = n1(0);
                coeff(0,1) = n1(1);
                coeff(0,2) = n1(2);
                coeff(1,0) = n2(0);
                coeff(1,1) = n2(1);
                coeff(1,2) = n2(2);
                coeff(2,0) = s(0);
                coeff(2,1) = s(1);
                coeff(2,2) = s(2);

                // Creo il vettore dei termini noti
                Vector3d term;
                term(0) = d1;
                term(1) = d2;
                term(2) = 0;

                // Risolvo il sistema e ottengo il punto Point
                Vector3d point;
                point = coeff.inverse() * term;

                // Vec3d Point;
                // Point.x = point[0];
                // Point.y = point[1];
                // Point.z = point[2];

                // Dal vettore s e dal punto Point posso ricavare la forma parametrica della retta



                // findTraces(s, Point, fracture, Id1, Id2);

                // Ciclo su tutti i vertici della prima frattura
                for(unsigned int n = 0; n < fracture.NumVertices[Id1]; n++)
                {
                    // Prima tratto il caso generale in cui collego un vertice con il successivo
                    if(n != fracture.NumVertices[Id1] - 1)
                    {
                        Vector3d v;
                        v(0) = fracture.vertices[Id1](0,n) - fracture.vertices[Id1](0,n+1);
                        v(1) = fracture.vertices[Id1](1,n) - fracture.vertices[Id1](1,n+1);
                        v(2) = fracture.vertices[Id1](2,n) - fracture.vertices[Id1](2,n+1);

                        Vector3d vec = s.cross(v);

                        Matrix3d A;
                        A << s, -v, vec;

                        Vector3d b;
                        b(0) = fracture.vertices[Id1](0,n) - point(0);
                        b(1) = fracture.vertices[Id1](1,n) - point(1);
                        b(2) = fracture.vertices[Id1](2,n) - point(2);

                        Vector3d sol;
                        sol = A.inverse() * b;

                        double t = sol(0);
                        double u = sol(1);

                        Vector3d intersection = point + t * s;

                        Vector3d point1;
                        point1 << fracture.vertices[Id1](0,n), fracture.vertices[Id1](1,n), fracture.vertices[Id1](2,n);
                        Vector3d verify = point1 + u * v;

                        double tol = 10e-16;

                        if((intersection - verify).norm() < tol)
                        {
                            cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                            cout << "Vertices: " << n << ", " << n+1 << endl;
                            cout << "Intersection: " << endl;
                            for(unsigned int i = 0; i < 3; i++)
                                cout << intersection(i) << endl;
                        }
                        else
                        {
                            cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                            cout << "Vertices: " << n << ", " << n+1 << endl;
                            cout << "Non si intersecano" << endl;
                        }
                    }

                    // Poi tratto il caso particolare del segmento formato dall'ultimo vertice con il primo
                    else
                    {
                        Vector3d v;
                        v(0) = fracture.vertices[Id1](0,n) - fracture.vertices[Id1](0,0);
                        v(1) = fracture.vertices[Id1](1,n) - fracture.vertices[Id1](1,0);
                        v(2) = fracture.vertices[Id1](2,n) - fracture.vertices[Id1](2,0);

                        Vector3d vec = s.cross(v);

                        Matrix3d A;
                        A << s, -v, vec;

                        Vector3d b;
                        b(0) = fracture.vertices[Id1](0,n) - point(0);
                        b(1) = fracture.vertices[Id1](1,n) - point(1);
                        b(2) = fracture.vertices[Id1](2,n) - point(2);

                        Vector3d sol;
                        sol = A.inverse() * b;

                        double t = sol(0);
                        double u = sol(1);

                        Vector3d intersection = point + t * s;

                        Vector3d point1;
                        point1 << fracture.vertices[Id1](0,n), fracture.vertices[Id1](1,n), fracture.vertices[Id1](2,n);
                        Vector3d verify = point1 + u * v;

                        double tol = 10e-16;

                        if((intersection - verify).norm() < tol)
                        {
                            cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                            cout << "Vertices: " << n << ", " << 0 << endl;
                            cout << "Intersection: " << endl;
                            for(unsigned int i = 0; i < 3; i++)
                                cout << intersection(i) << endl;
                        }
                        else
                        {
                            cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                            cout << "Vertices: " << n << ", " << 0 << endl;
                            cout << "Non si intersecano" << endl;
                        }
                    }
                }
            }
        }
    }
    return true;
}


}
