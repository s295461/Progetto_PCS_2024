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

// Funzione che trova le tracce
//bool findTraces(const Vec3d s, const Vec3d Point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2)
//{

//}

BoundingBox BBox3D(const MatrixXd& vertices) {
    BoundingBox bbox;
    //Inizializzato a valori infinitamente grandi positivi e negativi
    bbox.min = Vector3d::Constant(numeric_limits<double>::infinity());
    bbox.max = Vector3d::Constant(-numeric_limits<double>::infinity());
    //Aggiorna quando trova un nuovo massimo o minimo elemento per elemento
    for (int i = 0; i < vertices.cols(); ++i) {
        bbox.min = bbox.min.cwiseMin(vertices.col(i));
        bbox.max = bbox.max.cwiseMax(vertices.col(i));
    }

    return bbox;
}
// Questa funzione calcola l'intersezione tra due piani creati da due fratture
bool PlaneIntersection(const DiscreteFractureNetwork fracture)
{
    vector<FractureBBox> BBoxVect(fracture.numFracture);

    for(unsigned int i = 0; i < fracture.numFracture; i++){
        unsigned int ID = fracture.fractureID[i];
        BBoxVect[i].bbox = BBox3D(fracture.vertices[ID]);
        BBoxVect[i].fractureID = ID;
    }

    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        for(unsigned int j = 0; j < fracture.numFracture; j++)
        {
            if(i != j)
            {
                //Bisogna controllare se le bounding box si sovrappongono
                bool intersezione = (BBoxVect[i].bbox.min.x()> BBoxVect[j].bbox.max.x()) ||
                                    (BBoxVect[i].bbox.min.y()> BBoxVect[j].bbox.max.y()) ||
                                    (BBoxVect[i].bbox.min.z()> BBoxVect[j].bbox.max.z()) ||
                                    (BBoxVect[i].bbox.min.x()> BBoxVect[j].bbox.max.x()) ||
                                    (BBoxVect[i].bbox.min.y()> BBoxVect[j].bbox.max.y()) ||
                                    (BBoxVect[i].bbox.min.z()> BBoxVect[j].bbox.max.z());
                if (!intersezione){
                    continue;
                }
                // Seleziono due fratture
                unsigned int Id1 = BBoxVect[i].fractureID;
                unsigned int Id2 = BBoxVect[j].fractureID;


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

                // Ciclo su tutti i vertici della prima frattura

                for(unsigned int n = 0; n < fracture.NumVertices[Id1]; n++)
                {
                    unsigned int k = (n+1) % fracture.NumVertices[Id1];
                    Vector3d v;
                    v = fracture.vertices[Id1].col(k) - fracture.vertices[Id1].col(n);


                    for(unsigned int m = 0; m < fracture.NumVertices[Id2]; ++m){
                        unsigned int l = (m+1) % fracture.NumVertices[Id2];
                        Vector3d w,z;
                        w = fracture.vertices[Id2].col(l) - fracture.vertices[Id2].col(m);
                        z = fracture.vertices[Id2].col(m) - fracture.vertices[Id1].col(n);
                        //Coefficienti dell'equazione parametrica
                        double a = v.dot(v);
                        double b = v.dot(w);
                        double c = w.dot(w);
                        double d = v.dot(z);
                        double e = w.dot(z);

                        //Denominatore per il sistema di equazioni
                        double den = a*c - b*b;

                        double tol = 1e-16;
                        // Verifico se i segmenti sono paralleli
                        if (abs(den)< tol){
                            continue;
                        }

                        double t = (b*e-c*d)/den;
                        double u = (a*e-b*d)/den;

                        if (t>=0 && t<=1 && u>= 0 && u<=1){
                               Vector3d intersection = fracture.vertices[Id1].col(n) + t*v;

                                cout << "Intersezione trovata tra i poligoni:" << endl;
                                cout << "Frattura Id1: " << BBoxVect[i].fractureID << ", vertici [" << n << ", " << k << "]" << endl;
                                cout << "Frattura Id2: " << BBoxVect[j].fractureID <<  ", vertici [" << m << ", " << l << "]" << endl;
                                cout << "Punto di intersezione: " << intersection.transpose() << endl;
                            }

                    }
                }
            }
        }
    }
    return true;
}


}
