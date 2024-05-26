#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

namespace FractureNetwork {


bool ImportFracture(const string& filePath, DiscreteFractureNetwork& fracture, Traces trace)
{
    cout << "FR3_data" << endl;
    string fileNameFR3 = "/FR3_data.txt";
    if(!readFracture(filePath, fileNameFR3, fracture))
        return false;


    PlaneIntersection(fracture, trace);

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
        verticesFracture << x.transpose(), y.transpose(), z.transpose();
        //for(unsigned int n = 0; n < fracture.NumVertices[i]; n++)
        //{
        //    verticesFracture(0, n) = x[n];
        //    verticesFracture(1, n) = y[n];
        //    verticesFracture(2, n) = z[n];
        //}

        // Salvo la matrice di coordinate in un vettore di matrici
        fracture.vertices.push_back(verticesFracture);
    }

    // Chiudo il file__
    file.close();
    return true;
}


// Questa funzione calcola l'intersezione tra due piani creati da due fratture
bool PlaneIntersection(const DiscreteFractureNetwork fracture, Traces trace)
{

    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        for(unsigned int j = 0; j < fracture.numFracture; j++)
        {
            if(i < j)
            {
                // Seleziono due fratture
                unsigned int Id1 = fracture.fractureID[i];
                unsigned int Id2 = fracture.fractureID[j];

                /// BOUNDING BOX
                // if(bounding box = True)
                // faccio tutto il resto che c'è sotto, altimenti avanzo nel ciclo

                // Calcolo i vettori
                Vector3d u1 = fracture.vertices[Id1].col(0) - fracture.vertices[Id1].col(2);
                Vector3d v1 = fracture.vertices[Id1].col(1) - fracture.vertices[Id1].col(2);
                Vector3d u2 = fracture.vertices[Id2].col(0) - fracture.vertices[Id2].col(2);
                Vector3d v2 = fracture.vertices[Id2].col(1) - fracture.vertices[Id2].col(2);

                // Calcolo le normali ai due piani
                Vector3d n1 = u1.cross(v1) / (u1.norm() * v1.norm());
                Vector3d n2 = u2.cross(v2) / (u2.norm() * v2.norm());

                // Calcolo il vettore direzione della retta di intersezione
                Vector3d s = n1.cross(n2);

                // Calcolo il termine noto
                double d1;
                d1 = n1(0) * fracture.vertices[Id1](0,2) + n1(1) * fracture.vertices[Id1](1,2) + n1(2) * fracture.vertices[Id1](2,2);

                double d2;
                d2 = n2(0) * fracture.vertices[Id2](0,2) + n2(1) * fracture.vertices[Id2](1,2) + n2(2) * fracture.vertices[Id2](2,2);

                // Inserisco i coefficienti del sistema formato dai due piani ricavati dalle fratture e dal piano con normale s in una matrice
                Matrix3d coeff;
                coeff << n1, n2, s;

                // Creo il vettore dei termini noti
                Vector3d term;
                term(0) = d1;
                term(1) = d2;
                term(2) = 0;

                // Risolvo il sistema e ottengo il punto Point
                Vector3d point;
                point = coeff.transpose().inverse() * term;

                // Dal vettore s e dal punto Point posso ricavare la forma parametrica della retta



                findTraces(s, point, fracture, Id1, Id2, trace);

                // Ciclo su tutti i vertici della prima frattura
                // for(unsigned int n = 0; n < fracture.NumVertices[Id1]; n++)
                // {
                //     // Ciclo su tutti i vertici della seconda frattura
                //     for(unsigned int m = 0; m < fracture.NumVertices[Id2]; m++)
                //     {
                //         unsigned int k = (n+1) % fracture.NumVertices[Id1];
                //         unsigned int h = (m+1) % fracture.NumVertices[Id2];

                //         Vector3d v1;
                //         v1 = fracture.vertices[Id1].col(n) - fracture.vertices[Id1].col(k);
                //         Vector3d P1;
                //         P1 = fracture.vertices[Id1].col(k);

                //         Vector3d v2;
                //         v2 = fracture.vertices[Id2].col(m) - fracture.vertices[Id2].col(h);
                //         Vector3d P2;
                //         P2 = fracture.vertices[Id2].col(h);

                //         // Verifico che le rette non siano parallele alla retta di intersezione tra i piani
                //         if((v1.cross(s)).norm() != 0 && (v2.cross(s)).norm() != 0)
                //         {
                //             double t1;
                //             double u1;
                //             t1 = ((P1.cross(v1) - point.cross(v1)).dot(s.cross(v1))) / (((s.cross(v1)).norm())*((s.cross(v1)).norm()));
                //             u1 = ((point.cross(s) - P1.cross(s)).dot(v1.cross(s))) / (((v1.cross(s)).norm())*((v1.cross(s)).norm()));

                //             double t2;
                //             double u2;
                //             t2 = ((P2.cross(v2) - point.cross(v2)).dot(s.cross(v2))) / (((s.cross(v2)).norm())*((s.cross(v2)).norm()));
                //             u2 = ((point.cross(s) - P2.cross(s)).dot(v2.cross(s))) / (((v2.cross(s)).norm())*((v2.cross(s)).norm()));

                //             Vector3d intersection1 = P1 + u1 * v1;
                //             Vector3d verify1 = point + t1 * s;

                //             Vector3d intersection2 = P2 + u2 * v2;
                //             Vector3d verify2 = point + t2 * s;

                //             double tol = 1e-16;

                //             if((intersection1 - verify1).norm() < tol && u1 >= 0 && u1 <= 1)
                //             {
                //                 cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                //                 cout << "Vertices: " << n << ", " << k << endl;
                //                 cout << "Intersection: " << endl;
                //                 for(unsigned int i = 0; i < 3; i++)
                //                     cout << intersection1(i) << endl;
                //             }
                //             else
                //             {
                //                 cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                //                 cout << "Vertices: " << n << ", " << k << endl;
                //                 cout << "Non si intersecano" << endl;
                //             }

                //             if((intersection2 - verify2).norm() < tol && u2 >= 0 && u2 <= 1)
                //             {
                //                 cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                //                 cout << "Vertices: " << n << ", " << k << endl;
                //                 cout << "Intersection: " << endl;
                //                 for(unsigned int i = 0; i < 3; i++)
                //                     cout << intersection2(i) << endl;
                //             }
                //             else
                //             {
                //                 cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                //                 cout << "Vertices: " << n << ", " << k << endl;
                //                 cout << "Non si intersecano" << endl;
                //             }

                //         }
                //         else
                //         {
                //             cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                //             cout << "Vertices: " << n << ", " << k << endl;
                //             cout << "Sono parallele, non si intersecano" << endl;
                //         }

                //     }





                        // MatrixXd A(3,2);
                        // A(0,0) = s(0);
                        // A(1,0) = s(1);
                        // A(2,0) = s(2);
                        // A(0,1) = v(0);
                        // A(1,1) = v(1);
                        // A(2,1) = v(2);

                        // Vector3d b;
                        // b(0) = fracture.vertices[Id1](0,k) - point(0);
                        // b(1) = fracture.vertices[Id1](1,k) - point(1);
                        // b(2) = fracture.vertices[Id1](2,k) - point(2);

                        // Vector2d sol;
                        // // // sol = A.inverse() * b;
                        // sol = A.partialPivLu().solve(b);

                        // double t = sol(0);
                        // double u = sol(1);





                    //Vector3d vec = s.cross(v);

                    // MatrixXd A(3,2);
                    // A(0,0) = s(0);
                    // A(1,0) = s(1);
                    // A(2,0) = s(2);
                    // A(0,1) = v(0);
                    // A(1,1) = v(1);
                    // A(2,1) = v(2);

                    // Vector3d b;
                    // b(0) = fracture.vertices[Id1](0,k) - point(0);
                    // b(1) = fracture.vertices[Id1](1,k) - point(1);
                    // b(2) = fracture.vertices[Id1](2,k) - point(2);

                    // Vector2d sol;
                    // // sol = A.inverse() * b;
                    // sol = A.fullPivLu().solve(b);

                    // double t = sol(0);
                    // double u = sol(1);

                    // Vector3d point1;
                    // point1 << fracture.vertices[Id1](0,n), fracture.vertices[Id1](1,n), fracture.vertices[Id1](2,n);
                }
            }
        }
    // }
    return true;
}



// Funzione che trova le tracce
bool findTraces(const Vector3d s, const Vector3d point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2, Traces trace)
{
    unsigned int n = 0;
    unsigned int m = 0;

    vector<Vector3d> Point1;
    vector<Vector3d> Point2;


    double tol = 1e-10;

    // Nel ciclo while prendo ogni volta un segmento della prima e un segmento della seconda frattura
    while(n < fracture.NumVertices[Id1] && m < fracture.NumVertices[Id2])
    {
        // Calcolo gli indici del vertice successivo al vertice considerato, se considero l'ultimo vertice il successivo sarà lo 0.
        unsigned int k = (n+1) % fracture.NumVertices[Id1];
        unsigned int h = (m+1) % fracture.NumVertices[Id2];

        Vector3d v1;
        v1 = fracture.vertices[Id1].col(n) - fracture.vertices[Id1].col(k);
        Vector3d P1;
        P1 = fracture.vertices[Id1].col(k);

        Vector3d v2;
        v2 = fracture.vertices[Id2].col(m) - fracture.vertices[Id2].col(h);
        Vector3d P2;
        P2 = fracture.vertices[Id2].col(h);

        // Verifico che le rette non siano parallele alla retta di intersezione tra i piani
        if((v1.cross(s)).norm() != 0 || (v2.cross(s)).norm() != 0)
        {
            double t1;
            double u1;
            t1 = ((P1.cross(v1) - point.cross(v1)).dot(s.cross(v1))) / (((s.cross(v1)).norm())*((s.cross(v1)).norm()));
            u1 = ((point.cross(s) - P1.cross(s)).dot(v1.cross(s))) / (((v1.cross(s)).norm())*((v1.cross(s)).norm()));

            double t2;
            double u2;
            t2 = ((P2.cross(v2) - point.cross(v2)).dot(s.cross(v2))) / (((s.cross(v2)).norm())*((s.cross(v2)).norm()));
            u2 = ((point.cross(s) - P2.cross(s)).dot(v2.cross(s))) / (((v2.cross(s)).norm())*((v2.cross(s)).norm()));

            Vector3d intersection1 = point + t1 * s;
            Vector3d verify1 = P1 + u1 * v1;

            Vector3d intersection2 = point + t2 * s;
            Vector3d verify2 = P2 + u2 * v2;

            if((intersection1 - verify1).norm() < tol && u1 >= 0 && u1 <= 1)
            {
                cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                cout << "Vertices_1: " << n << ", " << k << endl;
                cout << "Intersection_1: " << endl;
                for(unsigned int i = 0; i < 3; i++)
                cout << intersection1(i) << endl;

                Point1.push_back(intersection1);
            }
            else
            {
                cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                cout << "Vertices_1: " << n << ", " << k << endl;
                cout << "Non si intersecano" << endl;
            }

            if((intersection2 - verify2).norm() < tol && u2 >= 0 && u2 <= 1)
            {
                cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                cout << "Vertices_2: " << m << ", " << h << endl;
                cout << "Intersection_2: " << endl;
                for(unsigned int i = 0; i < 3; i++)
                    cout << intersection2(i) << endl;

                Point2.push_back(intersection2);
            }
            else
            {
                cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                cout << "Vertices_2: " << m << ", " << h << endl;
                cout << "Non si intersecano" << endl;
            }
        }
        else
        {
            cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
            cout << "Vertices: " << n << ", " << k << endl;
            cout << "Sono parallele, non si intersecano" << endl;
        }

        n++;
        m++;
    }


    unsigned int idTrace = 0;

    // Se i punti coincidono ho già trovato la traccia e la salvo
    // CRASH PERCHE' ALCUNE INTERSEZIONI NON DANNO UN PUNTO COME RISULTATO, DUNQUE DA ERRORE LA DIFFERENZA NEL IF, DOVREBBE ANDARE APOSTO CON IL BOUNDING BOX
    if(((Point1[0] - Point2[0]).norm() < tol && (Point1[1] - Point2[1]).norm() < tol) || ((Point1[0] - Point2[1]).norm() < tol && (Point1[1] - Point2[0]).norm() < tol))
    {
        cout << "Coincidono" << endl;
        trace.traceId.push_back(idTrace);
        Vector2i id;
        id << Id1, Id2;
        trace.fractureId.push_back(id);
        MatrixXd verticesTrace(3,2);
        verticesTrace << Point1[0], Point1[1];
        trace.coordinates.push_back(verticesTrace);
        trace.numTraces++;

        idTrace++;
    }

    else
    {


        cout << "non coincidono" << endl;
    }


return true;
}






}
