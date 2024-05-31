#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

namespace FractureNetwork {


bool ImportFracture(const string& filePathInput, const string filePathOutput, DiscreteFractureNetwork& fracture, Traces& trace)
{
    // cout << "FR3_data" << endl;
    string fileNameFR3 = "/FR3_data.txt";
    if(!ReadFracture(filePathInput, fileNameFR3, fracture))
    {
        cerr << "Something wrong with the reading of the fracture" << endl;
        return false;
    }

    if(!FractureIntersection(fracture, trace))
    {
        cerr << "Something wrong while computing fracture intersection" << endl;
        return false;
    }

    string outputFileNameFR3 = "/FR3_traces.txt";
    if(!PrintOnFile(outputFileNameFR3, filePathOutput, trace))
    {
        cerr << "Something wrong while printing the result in a file" << endl;
        return false;
    }

    if(!TraceReorder(fracture, trace, filePathOutput))
    {
        cerr << "Something wrong while reordering traces" << endl;
        return false;
    }

    string reorderedOutputFileNameFR3 = "/FR3_traces_reordered.txt";
    if(!printTraces(reorderedOutputFileNameFR3, filePathOutput, trace, fracture))
    {
        cerr << "Something vrong printing the traces reordered" << endl;
        return false;
    }

    //cout << endl;
    //cout << "FR10_data" << endl;
    //string fileNameFR10 = "/FR10_data.txt";
    //readFracture(filePath, fileNameFR10, fracture);
    return true;
}


// Questa funzione apre un file, legge il contenuto e lo salva in strutture dati adeguate
bool ReadFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture)
{
    // Apro il file
    ifstream file;
    file.open(filePath + fileName);

    // Verifico che si apra correttamente
    if(file.fail())
    {
        cerr << "Error opening the file" << endl;
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
        // cout << "id: " << id << endl;
        // cout << "num: " << num << endl;

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

        // // Stampo i vertici
        // for(unsigned int n = 0; n < fracture.NumVertices[i]; n++)
        // {
        //     cout << scientific << setprecision(16) << "Vertices " << n+1 << ": x: " << x[n] << "; y: " << y[n] << "; z: " << z[n] << endl;
        // }

        // Salvo in una matrice le coordinate dei vertici della frattura
        MatrixXd verticesFracture(3,fracture.NumVertices[i]);
        verticesFracture << x.transpose(), y.transpose(), z.transpose();

        // Salvo la matrice di coordinate in un vettore di matrici
        fracture.vertices.push_back(verticesFracture);
    }

    // Chiudo il file__
    file.close();
    return true;
}


BoundingBox BBox3D(const MatrixXd& vertices)
{
    BoundingBox bbox;
    //Inizializzato a valori infinitamente grandi positivi e negativi
    bbox.min = Vector3d::Constant(-numeric_limits<double>::infinity());
    bbox.max = Vector3d::Constant(numeric_limits<double>::infinity());
    //Aggiorna quando trova un nuovo massimo o minimo elemento per elemento
    for (int i = 0; i < vertices.cols(); ++i) {
        bbox.min = bbox.min.cwiseMin(vertices.col(i));
        bbox.max = bbox.max.cwiseMax(vertices.col(i));
    }


    return bbox;
}


// Questa funzione calcola l'intersezione tra due piani creati da due fratture
bool FractureIntersection(const DiscreteFractureNetwork fracture, Traces& trace)
{
    vector<FractureBBox> BBoxVect(fracture.numFracture);

    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        unsigned int ID = fracture.fractureID[i];
        BBoxVect[i].bbox = BBox3D(fracture.vertices[ID]);
        BBoxVect[i].fractureID = ID;
    }
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        for(unsigned int j = 0; j < fracture.numFracture; j++)
        {
            if(i < j)
            {
                bool intersezione = !(BBoxVect[i].bbox.min.x() > BBoxVect[j].bbox.max.x() ||
                                      BBoxVect[i].bbox.min.y() > BBoxVect[j].bbox.max.y() ||
                                      BBoxVect[i].bbox.min.z() > BBoxVect[j].bbox.max.z() ||
                                      BBoxVect[i].bbox.max.x() < BBoxVect[j].bbox.min.x() ||
                                      BBoxVect[i].bbox.max.y() < BBoxVect[j].bbox.min.y() ||
                                      BBoxVect[i].bbox.max.z() < BBoxVect[j].bbox.min.z());
                if (!intersezione) {
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

                // Chiamo ora la funzione findTraces che trova la traccia tra la frattura con Id1 e la frattura con Id2.
                if(!FindTraces(s, point, fracture, Id1, Id2, trace))
                {
                    cerr << "Something wrong while searching traces" << endl;
                    return false;
                }
            }
        }
    }
    return true;
}


// Funzione che trova le tracce
bool FindTraces(const Vector3d s, const Vector3d point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2, Traces& trace)
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

            // Se il punto di intersezione coincide con il punto di verifica, allora è corretto
            if((intersection1 - verify1).norm() < tol && u1 >= 0 && u1 <= 1)
            {
                // cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                // cout << "Vertices_1: " << n << ", " << k << endl;
                // cout << "Intersection_1: " << endl;
                // for(unsigned int i = 0; i < 3; i++)
                //     cout << intersection1(i) << endl;



                Point1.push_back(intersection1);
            }
            // Se i punti non coincidono non si intersecano
            // else
            // {
            //     cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
            //     cout << "Vertices_1: " << n << ", " << k << endl;
            //     cout << "Non si intersecano" << endl;
            // }

            if((intersection2 - verify2).norm() < tol && u2 >= 0 && u2 <= 1)
            {
                // cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
                // cout << "Vertices_2: " << m << ", " << h << endl;
                // cout << "Intersection_2: " << endl;
                // for(unsigned int i = 0; i < 3; i++)
                //     cout << intersection2(i) << endl;

                Point2.push_back(intersection2);
            }
            // else
            // {
            //     cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
            //     cout << "Vertices_2: " << m << ", " << h << endl;
            //     cout << "Non si intersecano" << endl;
            // }
        }
        // else
        // {
        //     cout << "Id1: " << Id1 << ", Id2: " << Id2 << endl;
        //     cout << "Vertices: " << n << ", " << k << endl;
        //     cout << "Sono parallele, non si intersecano" << endl;
        // }

        // Passo ai due segmenti successivi
        n++;
        m++;
    }

    // Dopo aver girato su tutti i lati delle fratture ottengo un vettore di 4 punti, due saranno dati dall'intersezione della retta risultante dall'intersezione
    // tra i piani con la frattura Id1 e gli altri due dall'intersezione della retta con la frattura Id2.
    if(!Point1.empty() && !Point2.empty())
    {
        vector<Vector3d> intersectionPoints = {Point1[0], Point1[1], Point2[0], Point2[1]};

        // Ottengo le posizioni relative dei quattro punti sulla retta di intersezione tra i due piani
        double a = ((intersectionPoints[0] - point).dot(s)) / (s.norm() * s.norm());
        double b = ((intersectionPoints[1] - point).dot(s)) / (s.norm() * s.norm());
        double c = ((intersectionPoints[2] - point).dot(s)) / (s.norm() * s.norm());
        double d = ((intersectionPoints[3] - point).dot(s)) / (s.norm() * s.norm());

        if(a > b)
        {
            double a1 = a;
            a = b;
            b = a1;
        }

        if(c > d)
        {
            double c1 = c;
            c = d;
            d = c1;
        }

        if(c > a && c < b && b < d)
        {
            SaveTraces(c, b, point, s, trace, Id1, Id2);
            trace.Tips.push_back(true);
        }
        else if(a > c && a < d && d < b)
        {
            SaveTraces(a, d, point, s, trace, Id1, Id2);
            trace.Tips.push_back(true);
        }
        // else if((b > a && b < c && c < d) || (d > c && d < a && a < b))
        // {
        //     cout << "Non si intersecano" << endl;
        // }
        else if(c > a && c < d && d < b)
        {
            SaveTraces(c, d, point, s, trace, Id1, Id2);
            trace.Tips.push_back(false);
        }
        else if(a > c && a < b && b < d)
        {
            SaveTraces(a, b, point, s, trace, Id1, Id2);
            trace.Tips.push_back(false);
        }
        else if(abs((a - c)) <= tol && (b - d) <= tol)
        {
            SaveTraces(a, b, point, s, trace, Id1, Id2);
            trace.Tips.push_back(false);
        }
        else if(abs((a - c)) <= tol && b < d)
        {
            SaveTraces(a, b, point, s, trace, Id1, Id2);
            trace.Tips.push_back(false);
        }
        else if(abs((a - c)) <= tol && d < b)
        {
            SaveTraces(c, d, point, s, trace, Id1, Id2);
            trace.Tips.push_back(false);
        }
        else if(abs((b - d)) <= tol && c < a)
        {
            SaveTraces(a, b, point, s, trace, Id1, Id2);
            trace.Tips.push_back(false);
        }
        else if(abs((b - d)) <= tol && a < c)
        {
            SaveTraces(c, d, point, s, trace, Id1, Id2);
            trace.Tips.push_back(false);
        }
    }
    return true;
}


// Con questa funzione salvo tutti i valori legati ad una traccia nella struttura Traces.
void SaveTraces(double n, double m, Vector3d point, Vector3d s, Traces& trace, unsigned int Id1, unsigned int Id2)
{
    MatrixXd Coordinates(3,2);
    Vector3d P,Q;
    P = point+s*n;
    Q = point+s*m;
    Coordinates << P,Q;
    trace.coordinates.push_back(Coordinates);
    Vector2i idFracture;
    idFracture << Id1, Id2;
    trace.fractureId.push_back(idFracture);
    trace.traceId.push_back(trace.numTraces);
    double lengthSegment = PointDistance(P, Q);
    trace.length.push_back(lengthSegment);
    trace.numTraces++;
}



// Questa funzione stampa i risultati in un file
bool PrintOnFile(const string fileName, const string filePath, Traces trace)
{
    ofstream file;
    file.open(filePath + fileName);
    if(file.fail())
    {
        cerr << "Error opening the file" << endl;
        return false;
    }

    file << "# Number of traces" << endl;
    file << trace.numTraces << endl;

    file << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for(unsigned int i = 0; i < trace.numTraces; i++)
    {
        file << scientific << setprecision(16) << trace.traceId[i] << "; " << trace.fractureId[i][0] << "; " << trace.fractureId[i][1] << "; " <<
            trace.coordinates[i](0,0) << "; " << trace.coordinates[i](1,0) << "; " << trace.coordinates[i](2,0) <<
            "; " << trace.coordinates[i](0,1) << "; " << trace.coordinates[i](1,1) << "; " << trace.coordinates[i](2,1) << endl;
    }
    file.close();
    return true;
}


bool TraceReorder(DiscreteFractureNetwork& fracture, Traces& trace, const string filePathOutput)
{
    double tol = 1e-10;
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        vector<unsigned int> passing = {};
        vector<unsigned int> notPassing = {};
        vector<double> lengthPassing = {};
        vector<double> lengthNotPassing = {};
        for(unsigned int j = 0; j < trace.numTraces; j++)
        {
            // Cerco nel vettore di vettori fractureId le posizioni in cui una delle due fratture che generano la traccia corrisponde alla frattura che sto considerando.
            if(find(trace.fractureId[j].data(), trace.fractureId[j].data() + trace.fractureId[j].size(), fracture.fractureID[i]) != trace.fractureId[j].data() + trace.fractureId[j].size())
            {
                unsigned int position = j;
                unsigned int counter = 0;
                unsigned int n = 0;
                // Ciclo su tutti i segmenti della frattura e verifico se uno dei vertici della traccia appartiene al segmento.
                while(n < fracture.NumVertices[i])
                {
                    unsigned int k = (n+1) % fracture.NumVertices[i];
                    // PQ è il vettore formato da due vertici consecutivi della frattura.
                    Vector3d PQ = fracture.vertices[i].col(n) - fracture.vertices[i].col(k);
                    // PA è il vettore formato da un vertice della frattura con un estremo della traccia.
                    Vector3d PA = fracture.vertices[i].col(n) - trace.coordinates[position].col(0);
                    // PB è il segmento formato dal vertice della frattura con l'altro vertice della traccia.
                    Vector3d PB = fracture.vertices[i].col(n) - trace.coordinates[position].col(1);
                    // Se il prodotto vettoriale tra PQ e PA oppure tra PQ e PB è zero, allora A o B appartengono al segmento PQ, dunque incremento di uno il contatore.
                    if((PQ.cross(PA)).norm() <= tol || (PQ.cross(PB)).norm() <= tol)
                        counter++;
                    n++;
                }
                // Finito il ciclo su tutti i segmenti, se il contatore è uguale a 2 vuol dire che i vertici della traccia appartengono a due segmenti della frattura,
                // dunque la traccia è passante per quella frattura. Memorizzo gli le posizioni, che corrispondono anche agli id delle tracce, in due liste, una per le tracce
                // passanti e una per quelle non passanti.
                if(counter == 2)
                    passing.push_back(position);
                else
                    notPassing.push_back(position);
            }
        }
        // Creo altre due liste, una per la lunghezza delle tracce passanti e una per la lunghezza delle tracce non passanti e inserisco le lunghezze mantenendo l'ordine degli id delle tracce.
        for(unsigned int a : passing)
            lengthPassing.push_back(trace.length[a]);
        for(unsigned int a : notPassing)
            lengthNotPassing.push_back(trace.length[a]);

        // Ora riordino la lista delle lunghezze in ordine decrescente scambiando anche le posizioni della lista di indici, così da mantenere la corrispondenza.
        if(!reordering(passing, lengthPassing))
        {
            cerr << "Error reordering the traces" << endl;
            return false;
        }

        if(!reordering(notPassing, lengthNotPassing))
        {
            cerr << "Error reordering the traces" << endl;
            return false;
        }
        trace.passingReordered.push_back(passing);
        trace.notPassingReordered.push_back(notPassing);


    }
    return true;
}


double PointDistance(Vector3d P, Vector3d Q)
{
    return sqrt(((P[0]-Q[0]) * (P[0] - Q[0])) + ((P[1] - Q[1]) * (P[1] - Q[1])) + ((P[2] - Q[2]) * (P[2] - Q[2])));
}


bool reordering(vector<unsigned int>& idTraces, vector<double>& length)
{
    if(idTraces.size() != length.size())
        return false;
    vector<pair<double, unsigned int>> couple;
    for(unsigned int i = 0; i < idTraces.size(); i++)
        couple.push_back((make_pair(length[i], idTraces[i])));

    sort(couple.begin(), couple.end(), [](const pair<double, unsigned int>& a, const pair<double, unsigned int>& b){
        return a.first > b.first;
    });
    for(size_t i = 0; i < couple.size(); i++)
    {
        length[i] = couple[i].first;
        idTraces[i] = couple[i].second;
    }
    return true;
}


// bool printTraces(const string fileName, const string filePath, Traces trace, DiscreteFractureNetwork fracture)
// {
//     ofstream file;
//     file.open(filePath + fileName);
//     if(file.fail())
//     {
//         cerr << "Error opening the file" << endl;
//         return false;
//     }

//     for(unsigned int i = 0; i < fracture.numFracture; i++)
//     {
//         file << "# FractureId; NumTraces" << endl;
//         file << fracture.fractureID[i] << "; " << trace.numTraces << endl;
//         file << "# TraceId; Tips; Length" << endl;
//         for(unsigned int j = 0; j < trace.passingReordered.size(); j++)
//             file << trace.passingReordered[i][j] << trace.Tips[trace.passingReordered[i][j]] << trace.length[trace.passingReordered[i][j]] << endl;
//         for(unsigned int j = 0; j < trace.notPassingReordered.size(); j++)
//             file << trace.notPassingReordered[i][j] << trace.Tips[trace.notPassingReordered[i][j]] << trace.length[trace.notPassingReordered[i][j]] << endl;
//     }

//     file.close();
//     return true;
// }




}
