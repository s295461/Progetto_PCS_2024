#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <string>

using namespace std;

namespace FractureNetwork {


bool ImportFracture(const string fileNameInput, const string fileNameOutput, const string fileNameOutputReordered,
                    const string filePathInput, const string filePathOutput , DiscreteFractureNetwork& fracture, Traces& trace)
{
    if(!ReadFracture(filePathInput, fileNameInput, fracture))
    {
        cerr << "Something wrong with the reading of the fracture" << endl;
        return false;
    }

    if(!FractureIntersection(fracture, trace))
    {
        cerr << "Something wrong while computing fracture intersection" << endl;
        return false;
    }

    if(!PrintOnFile(fileNameOutput, filePathOutput, trace))
    {
        cerr << "Something wrong while printing the result in a file" << endl;
        return false;
    }


    if(!TraceReorder(fracture, trace))
    {
        cerr << "Something wrong while reordering traces" << endl;
        return false;
    }


    if(!printTraces(fileNameOutputReordered, filePathOutput, trace, fracture))

    {
        cerr << "Something vrong printing the traces reordered" << endl;
        return false;
    }

    return true;
}



// Questa funzione svuota tutti gli elementi all'interno della struttura DiscreteFractureNetwork
void clearDiscreteFractureNetwork(DiscreteFractureNetwork& fracture)
{
    fracture.numFracture = 0;
    fracture.fractureID.clear();
    fracture.NumVertices.clear();
    fracture.vertices.clear();
}

// Questa funzione svuota tutti gli elementi all'interno della struttura Traces
void clearTraces(Traces& trace)
{
    trace.numTraces = 0;
    trace.traceId.clear();
    trace.fractureId.clear();
    trace.coordinates.clear();
    trace.length.clear();
    trace.traceReordered.clear();
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

                Point1.push_back(intersection1);

            if((intersection2 - verify2).norm() < tol && u2 >= 0 && u2 <= 1)
                Point2.push_back(intersection2);
        }

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
            SaveTraces(c, b, point, s, trace, Id1, Id2);

        else if(a > c && a < d && d < b)
            SaveTraces(a, d, point, s, trace, Id1, Id2);

        else if(c > a && c < d && d < b)
            SaveTraces(c, d, point, s, trace, Id1, Id2);

        else if(a > c && a < b && b < d)
            SaveTraces(a, b, point, s, trace, Id1, Id2);

        else if(abs((a - c)) <= tol && (b - d) <= tol)
            SaveTraces(a, b, point, s, trace, Id1, Id2);

        else if(abs((a - c)) <= tol && b < d)
            SaveTraces(a, b, point, s, trace, Id1, Id2);

        else if(abs((a - c)) <= tol && d < b)
            SaveTraces(c, d, point, s, trace, Id1, Id2);

        else if(abs((b - d)) <= tol && c < a)
            SaveTraces(a, b, point, s, trace, Id1, Id2);

        else if(abs((b - d)) <= tol && a < c)
            SaveTraces(c, d, point, s, trace, Id1, Id2);
    }
    return true;
}

  
// Con questa funzione salvo tutti i valori legati ad una traccia nella struttura Traces.
void SaveTraces(double n, double m, Vector3d point, Vector3d s, Traces& trace, unsigned int Id1, unsigned int Id2)
{
    MatrixXd Coordinates(3,2);
    Vector3d P,Q;
    P = point+s*m;
    Q = point+s*n;
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

double PointDistance(Vector3d P, Vector3d Q){
    return sqrt(((P[0]-Q[0])*(P[0]-Q[0])+((P[1]-Q[1])*(P[1]-Q[1]))+((P[2]-Q[2])*(P[2]-Q[2]))));
}


bool TraceReorder(DiscreteFractureNetwork& fracture, Traces& trace)

{
    double tol = 1e-10;
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        vector<tuple<unsigned int, bool, double>> fractureTraces;
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

        for(unsigned int n = 0; n < passing.size(); n++)
        {
            tuple<unsigned int, bool, double> triplets;
            triplets = make_tuple(passing[n], false, lengthPassing[n]);
            fractureTraces.push_back(triplets);
        }
        for(unsigned int n = 0; n < notPassing.size(); n++)
        {
            tuple<unsigned int, bool, double> triplets;
            triplets = make_tuple(notPassing[n], true, lengthNotPassing[n]);
            fractureTraces.push_back(triplets);
        }

        trace.traceReordered.push_back(fractureTraces);
    }
    return true;
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



bool printTraces(const string fileName, const string filePath, Traces trace, DiscreteFractureNetwork fracture)

{
    ofstream file;
    file.open(filePath + fileName);
    if(file.fail())
    {
        cerr << "Error opening the file" << endl;
        return false;
    }

    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        unsigned int fractureNumTraces = trace.traceReordered[i].size();
        file << "# FractureId; NumTraces" << endl;
        file << fracture.fractureID[i] << "; " << fractureNumTraces << endl;
        file << "# TraceId; Tips; Length" << endl;
        if(fractureNumTraces != 0)
            for(unsigned int j = 0; j < fractureNumTraces; j++)
                file << get<0>(trace.traceReordered[i][j]) << "; " << get<1>(trace.traceReordered[i][j]) << "; " << get<2>(trace.traceReordered[i][j]) << endl;
        file << endl;

    }
    file.close();
    return true;
}

}

// **********************************************************************************************************************************************


namespace PolygonalMesh {

// Questa funzione effettua il taglio delle fratture in sottopoligoni e salva i risultati in una mesh poligonale.
bool fractureCut(FractureNetwork::DiscreteFractureNetwork& fracture, FractureNetwork::Traces& trace, Cell0D& Cell0D, Cell1D& Cell1D, Cell2D& Cell2D)
{
    double tol = 1e-10;
    // Prendo in esame una frattura alla volta
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        unsigned int id = fracture.fractureID[i];
        vector<pair<vector<Vector3d>, vector<pair<vector<Vector3d>, unsigned int>>>> subFracture;
        vector<pair<vector<Vector3d>, vector<pair<vector<Vector3d>, unsigned int>>>> subFracture1;
        vector<Vector3d> subfractureVertices;
        vector<vector<Vector3d>> subfractureVertices1;
        vector<Vector3d> subfractureVerticesDx;
        vector<Vector3d> subfractureVerticesSx;
        vector<pair<vector<Vector3d>, unsigned int>> subTracesDx;
        vector<pair<vector<Vector3d>, unsigned int>> subTracesSx;
        // Creo il secondo elemento di subFracture
        vector<pair<vector<Vector3d>, unsigned int>> coupleTraces;
        for(unsigned int n = 0; n < trace.traceReordered[id].size(); n++)
        {
            tuple<unsigned int, bool, double> triplets = trace.traceReordered[id][n];
            unsigned int traceId = get<0>(triplets);
            vector<Vector3d> coord;
            Vector3d point = trace.coordinates[traceId].col(0);
            coord.push_back(point);
            Vector3d point1 = trace.coordinates[traceId].col(1);
            coord.push_back(point1);
            coupleTraces.push_back(make_pair(coord, traceId));
        }
        // Prendo le coordinate della frattura iniziale
        for(unsigned int n = 0; n < fracture.vertices[id].cols(); n++)
        {
            Vector3d point = fracture.vertices[id].col(n);
            subfractureVertices.push_back(point);
        }
        // Creo subFracture iniziale
        subFracture.push_back(make_pair(subfractureVertices, coupleTraces));
        // Prendo una traccia alla volta
        for(unsigned int j = 0; j < trace.traceReordered[id].size(); j++)
        {
            // Prendo l'id della traccia che sto studiando
            tuple<unsigned int, bool, double> triplets = trace.traceReordered[id][j];
            unsigned int traceId = get<0>(triplets);
            subFracture1 = subFracture;
            // Ciclo su tutti gli elementi di subFracture, quindi tutte le sottofratture che ho già creato
            for(unsigned int n = 0; n < subFracture.size(); n++)
            {
                // Prendo le coordinate della sottofrattura da tagliare
                vector<Vector3d> fractureCoord = get<0>(subFracture[n]);
                // Verifico che la sottofrattura che sto considerando ha dei vertici
                if(fractureCoord.size() == 0)
                    continue;
                // Calcolo la normale alla frattura, mi servirà dopo
                Vector3d normal = (fractureCoord[1] - fractureCoord[0]).cross(fractureCoord[2] - fractureCoord[1]);
                // Prendo le informazioni sulle tracce che tagliano questa sottofrattura
                vector<pair<vector<Vector3d>,unsigned int>> fractureTraces = get<1>(subFracture[n]);
                vector<Vector3d> traceCoord;
                // Se la traccia che sto studiando è presente tra le tracce che tagliano la sottofrattura ne prendo le coordinate
                unsigned int m = 0;
                unsigned int flag = 1;
                while(m < fractureTraces.size())
                {
                    if(get<1>(fractureTraces[m]) == traceId)
                    {
                        traceCoord = get<0>(fractureTraces[m]);
                        flag = 0;
                        m = fractureTraces.size();
                    }
                    m++;
                }
                // Se l'id della traccia che sto studiando (j) non è presente tra le tracce che tagliano la sottofrattura (n) passo alla sottofrattura successiva.
                if(flag == 1)
                    continue;

                // Se ho superato questo if so che la sottofrattura che sto considerando verrà tagliata, dunque elimino la frattura iniziale da subFracture1
                auto it = remove(subFracture1.begin(), subFracture1.end(), subFracture[n]);
                subFracture1.erase(it, subFracture1.end());

                // Devo controllare se la traccia è passante o meno, se è non passante devo estenderla finchè non incontro i lati della sottofrattura
                if(get<1>(triplets))
                    traceCoord = extendTraces(fractureCoord, traceCoord);

                if(!createSubfracture(fractureCoord, traceCoord, subfractureVertices1))
                {
                    cerr << "Error creating subtraces" << endl;
                    return false;
                }
                Vector3d cuttingTrace = traceCoord[1] - traceCoord[0];
                // Divido le sottofratture in destra e sinistra
                for(unsigned int h = 0; h < subfractureVertices1.size(); h++)
                {
                    Vector3d vec;
                    // Giro sui vertici della sottofrattura h
                    for(unsigned int k = 0; k < subfractureVertices1[h].size(); k++)
                    {
                        // Entro nell'if solo se il punto che sto considerando, cioè quello con indice k, non è uno dei due estremi della traccia per cui ho tagliato.
                        if((subfractureVertices1[h][k] - traceCoord[0]).norm() > tol && (subfractureVertices1[h][k] - traceCoord[1]).norm() > tol)
                        {
                            // Memorizzo il vettore ed esco dal ciclo
                            vec = subfractureVertices1[h][k] - traceCoord[0];
                            k = subfractureVertices1[h].size();
                        }
                    }

                    // Calcolo il prodotto vettoriale tra la traccia che ho tagliato e il vettore tra un estremo della traccia e un vertice della frattura
                    Vector3d crossProduct = cuttingTrace.cross(vec);
                    // Se il prodotto scalare tra la normale alla frattura e il prodotto vettoriale è positivo, vuol dire che la normale alla frattura e il prodotto vettoriale hanno lo stesso verso, dunque la sottofrattura sta a sinistra della traccia che ha tagliato
                    if(normal.dot(crossProduct) > 0)
                        subfractureVerticesSx = subfractureVertices1[h];
                    // Altrimenti i due vettori hanno verso opposto, dunque la sottofrattura sta a destra della traccia che ha tagliato
                    else
                        subfractureVerticesDx = subfractureVertices1[h];
                }
                // Una volta che ho diviso i suoi elementi in dx e sx, svuoto subfractureVertices1
                subfractureVertices1.clear();

                // Per ogni traccia all'interno della sottofrattura
                for(unsigned int s = 0; s < fractureTraces.size(); s++)
                {
                    unsigned int traceId1 = get<1>(fractureTraces[s]);
                    // La considero solo se è diversa dalla traccia per cui ho tagliato
                    if(traceId1 != traceId)
                    {
                        // Prendo i vettori tra i due estremi della traccia che voglio smistare in destra o sinistra e un estremo della traccia per cui ho tagliato
                        Vector3d vec1 = get<0>(fractureTraces[s])[0] - traceCoord[0];
                        Vector3d vec2 = get<0>(fractureTraces[s])[1] - traceCoord[0];
                        // Faccio i prodotti vettoriali tra la traccia per cui ho tagliato e i due vettori
                        Vector3d crossProduct1 = cuttingTrace.cross(vec1);
                        Vector3d crossProduct2 = cuttingTrace.cross(vec2);
                        // Se i prodotti scalari tra la normale alla frattura e i due prodotti vettoriali sono entrambi positivi o nulli, significa che la traccia si trova a sinistra della traccia per cui ho tagliato
                        if(normal.dot(crossProduct1) >= 0 && normal.dot(crossProduct2) >= 0)
                            subTracesSx.push_back(fractureTraces[s]);
                        // Se i prodotti scalari tra la normale alla frattura e i due prodotti vettoriali sono entrambi negativi o nulli, significa che la traccia si trova a destra della traccia per cui ho tagliato
                        else if(normal.dot(crossProduct1) <= 0 && normal.dot(crossProduct2) <= 0)
                            subTracesDx.push_back(fractureTraces[s]);
                        // Se un prodotto scalare è positivo e l'altro è negativo vuol dire che le tracce si intersecano, trovo il punto di intersezione e divido le sottotracce in destra e sinistra
                        else if(normal.dot(crossProduct1) > 0 && normal.dot(crossProduct2) < 0)
                        {
                            Vector3d traceIntersection;
                            // Trovo il vettore della traccia che voglio smistare
                            Vector3d vectorTrace = get<0>(fractureTraces[s])[1] - get<0>(fractureTraces[s])[0];
                            // Calcolo il punto di intersezione
                            double t;
                            double u;
                            t = ((get<0>(fractureTraces[s])[0].cross(vectorTrace) - traceCoord[0].cross(vectorTrace)).dot(cuttingTrace.cross(vectorTrace))) / (((cuttingTrace.cross(vectorTrace)).norm())*((cuttingTrace.cross(vectorTrace)).norm()));
                            u = ((traceCoord[0].cross(cuttingTrace) - get<0>(fractureTraces[s])[0].cross(cuttingTrace)).dot(vectorTrace.cross(cuttingTrace))) / (((vectorTrace.cross(cuttingTrace)).norm())*((vectorTrace.cross(cuttingTrace)).norm()));
                            Vector3d intersection = traceCoord[0] + t * cuttingTrace;
                            Vector3d verify = get<0>(fractureTraces[s])[0] + u * vectorTrace;
                            // Se il punto di intersezione coincide con il punto di verifica, allora è corretto, lo salvo come punto di intersezione tra le due tracce
                            if((intersection - verify).norm() < tol && t >= 0 && t <= 1)
                                traceIntersection = intersection;
                            vector<Vector3d> newTraceDx = {traceIntersection, get<0>(fractureTraces[s])[1]};
                            pair<vector<Vector3d>,unsigned int> traceDx = make_pair(newTraceDx, traceId1);
                            subTracesDx.push_back(traceDx);
                            vector<Vector3d> newTraceSx = {traceIntersection, get<0>(fractureTraces[s])[0]};
                            pair<vector<Vector3d>,unsigned int> traceSx = make_pair(newTraceSx, traceId1);
                            subTracesSx.push_back(traceSx);
                        }
                        else if(normal.dot(crossProduct1) < 0 && normal.dot(crossProduct2) > 0)
                        {
                            Vector3d traceIntersection;
                            // Trovo il vettore della traccia che voglio smistare
                            Vector3d vectorTrace = get<0>(fractureTraces[s])[1] - get<0>(fractureTraces[s])[0];
                            // Calcolo il punto di intersezione
                            double t;
                            double u;
                            t = ((get<0>(fractureTraces[s])[0].cross(vectorTrace) - traceCoord[0].cross(vectorTrace)).dot(cuttingTrace.cross(vectorTrace))) / (((cuttingTrace.cross(vectorTrace)).norm())*((cuttingTrace.cross(vectorTrace)).norm()));
                            u = ((traceCoord[0].cross(cuttingTrace) - get<0>(fractureTraces[s])[0].cross(cuttingTrace)).dot(vectorTrace.cross(cuttingTrace))) / (((vectorTrace.cross(cuttingTrace)).norm())*((vectorTrace.cross(cuttingTrace)).norm()));
                            Vector3d intersection = traceCoord[0] + t * cuttingTrace;
                            Vector3d verify = get<0>(fractureTraces[s])[0] + u * vectorTrace;
                            // Se il punto di intersezione coincide con il punto di verifica, allora è corretto, lo salvo come punto di intersezione tra le due tracce
                            if((intersection - verify).norm() < tol && t >= 0 && t <= 1)
                                traceIntersection = intersection;
                            vector<Vector3d> newTraceDx = {traceIntersection, get<0>(fractureTraces[s])[0]};
                            pair<vector<Vector3d>,unsigned int> traceDx = make_pair(newTraceDx, traceId1);
                            subTracesDx.push_back(traceDx);
                            vector<Vector3d> newTraceSx = {traceIntersection, get<0>(fractureTraces[s])[1]};
                            pair<vector<Vector3d>,unsigned int> traceSx = make_pair(newTraceSx, traceId1);
                            subTracesSx.push_back(traceSx);
                        }
                    }
                }
                // Salvo le coppie sottofrattura-sottotracce in subFracture1
                subFracture1.push_back(make_pair(subfractureVerticesDx, subTracesDx));
                subFracture1.push_back(make_pair(subfractureVerticesSx, subTracesSx));
                // Svuoto le strutture che non mi servono più
                subfractureVerticesDx.clear();
                subfractureVerticesSx.clear();
                subTracesDx.clear();
                subTracesSx.clear();
            }
            // Una volta che ho esaminato tutte le coppie di subFracture, passo gli elementi di subFracture1 in subFracture e ricomincio il ciclo tagliando per la frattura successiva
            subFracture.clear();
            subFracture = subFracture1;
            subFracture1.clear();
        }
    }
    return true;
}



// Questa funzione effettua il taglio di una frattura per una traccia passante
bool createSubfracture(vector<Vector3d> subfracture, vector<Vector3d> cuttingTrace, vector<vector<Vector3d>>& subfractureVertices1)
{
    unsigned int numVertices = subfracture.size();
    double tol = 1e-10;
    vector<Vector3d> Points;
    unsigned int n = 0;
    unsigned int pos1 = 0;
    unsigned int pos2 = 0;
    // Creo un contatore che viene incrementato ogni volta cha aggiungo un elemento a Points
    unsigned int c = 0;
    // Con questo ciclo while creo un vettore Points che contiene i vertici della frattura e della traccia nell'ordine in cui li trovo percorrendo il perimetro della frattura.
    while(n < numVertices)
    {
        unsigned int k = (n + 1) % numVertices;
        Vector3d P = subfracture[n];
        Vector3d Q = subfracture[k];
        Vector3d A = cuttingTrace[0];
        Vector3d B = cuttingTrace[1];
        Vector3d PQ = Q - P;
        Vector3d PA = A - P;
        Vector3d PB = B - P;
        Points.push_back(P);
        c++;
        if((PQ.cross(PA)).norm() < tol)
        {
            Points.push_back(A);
            pos1 = c;
            c++;
        }
        else if((PQ.cross(PB)).norm() < tol)
        {
            Points.push_back(B);
            pos2 = c;
            c++;
        }
        n++;

    }

    vector<Vector3d> subfracture1, subfracture2;
    unsigned int m = 0;
    while(m < Points.size())
    {
        subfracture1.push_back(Points[m]);
        if(m == pos1)
        {
            m = pos2;
            subfracture1.push_back(Points[m]);
        }
        else if(m == pos2)
        {
            m = pos1;
            subfracture1.push_back(Points[m]);
        }
        m++;
    }
    unsigned int t = 0;
    if(pos1 < pos2)
    {
        t = pos1;
        while(t < Points.size())
        {
            subfracture2.push_back(Points[t]);
            if(t == pos2)
                t = Points.size() - 1;
            t++;
        }
    }
    else
    {
        t = pos2;
        while(t < Points.size())
        {
            subfracture2.push_back(Points[t]);
            if(t == pos1 )
                t = Points.size() - 1;
            t++;
        }
    }
    subfractureVertices1.push_back(subfracture1);
    subfractureVertices1.push_back(subfracture2);

    return true;
}

vector<Vector3d> extendTraces(vector<Vector3d> subFractureVertices, vector<Vector3d> subTraceVertices)
{
    double tol = 1e-10;
    Vector3d vec = subTraceVertices[1] - subTraceVertices[0];
    Vector3d point = subTraceVertices[0];
    vector<Vector3d> extendedVertices;
    unsigned int i = 0;
    while(i < subFractureVertices.size())
    {
        unsigned int k = (i + 1) % subFractureVertices.size();
        Vector3d vec1 = subFractureVertices[k] - subFractureVertices[i];
        Vector3d point1 = subFractureVertices[i];
        if((vec1.cross(vec)).norm() > tol)
        {
            double t;
            double u;
            t = ((point1.cross(vec1) - point.cross(vec1)).dot(vec.cross(vec1))) / (((vec.cross(vec1)).norm())*((vec.cross(vec1)).norm()));
            u = ((point.cross(vec) - point1.cross(vec)).dot(vec1.cross(vec))) / (((vec1.cross(vec)).norm())*((vec1.cross(vec)).norm()));


            Vector3d intersection = point + t * vec;
            Vector3d verify = point1 + u * vec1;

            // Se il punto di intersezione coincide con il punto di verifica, allora è corretto
            if((intersection - verify).norm() < tol && u >= 0 && u <= 1)
            {
                extendedVertices.push_back(intersection);
            }
        }
        i++;
    }
    return extendedVertices;

}

// // Questa funzione trova i punti di intersezione tra una traccia e tutte quelle che tagliano precedentemente la frattura
// vector<Vector3d> intersectTraces(FractureNetwork::Traces trace, unsigned int idTrace, unsigned int idFracture)
// {
//     // Memorizzo in un vettore traceBefore gli id delle tracce che tagliano la frattura prima di quella che sto studiando
//     vector<unsigned int> traceBefore;
//     unsigned int m = 0;
//     while(m < trace.traceReordered[idFracture].size())
//     {
//         tuple<unsigned int, bool, double> triplets = trace.traceReordered[idFracture][m];
//         unsigned int id = get<0>(triplets);
//         if(id == idTrace)
//             m = trace.traceReordered[idFracture].size();
//         traceBefore.push_back(id);
//         m++;
//     }
//     double tol = 1e-10;
//     vector<Vector3d> intersectionPoint;
//     unsigned int n = 1;
//     // Creo un vettore nullo
//     Vector3d null;
//     null << NULL, NULL, NULL;
//     // PROBLEMA: IO VOGLIO SOLO LE TRACCE DI QUELLA FRATTURA
//     while(n <= traceBefore.size())
//     {
//         Vector3d vec = trace.coordinates[idTrace].col(1) - trace.coordinates[idTrace].col(0);
//         Vector3d point = trace.coordinates[idTrace].col(0);
//         Vector3d vec1 = trace.coordinates[traceBefore[n]].col(1) - trace.coordinates[traceBefore[n]].col(0);
//         Vector3d point1 = trace.coordinates[traceBefore[n]].col(0);
//         if((vec1.cross(vec)).norm() != 0)
//         {
//             double t;
//             double u;
//             t = ((point1.cross(vec1) - point.cross(vec1)).dot(vec.cross(vec1))) / (((vec.cross(vec1)).norm())*((vec.cross(vec1)).norm()));
//             u = ((point.cross(vec) - point1.cross(vec)).dot(vec1.cross(vec))) / (((vec1.cross(vec)).norm())*((vec1.cross(vec)).norm()));

//             Vector3d intersection = point + t * vec;
//             Vector3d verify = point1 + u * vec1;

//             // Se il punto di intersezione coincide con il punto di verifica, allora è corretto, lo inserisco in un vettore nella posizione corrispondente alla traccia intersecata dalla traccia che sto esaminando
//             if((intersection - verify).norm() < tol)
//             {
//                 intersectionPoint.push_back(intersection);
//             }
//             // Se non si interseca inserisco un vettore nullo
//             else
//             {
//                 intersectionPoint.push_back(null);
//             }
//         }
//         n++;
//     }
//     return intersectionPoint;
// }

// // Questa funzione interseca tra di loro tutte le tracce di una frattura e trova i punti di intersezione
// vector<tuple<Vector3d, unsigned int, unsigned int>> intersectTraces(FractureNetwork::Traces trace, unsigned int idFracture)
// {
//     // Memorizzo in un vettore tracesId gli id delle tracce che tagliano la frattura che sto studiando
//     vector<unsigned int> tracesId;
//     for(unsigned int i = 0; i < trace.traceReordered[idFracture].size(); i++)
//     {
//         tuple<unsigned int, bool, double> triplets = trace.traceReordered[idFracture][i];
//         unsigned int id = get<0>(triplets);
//         tracesId.push_back(id);
//     }
//     double tol = 1e-10;
//     vector<tuple<Vector3d, unsigned int, unsigned int>> traceIntersection;
//     for(unsigned int n = 0; n < tracesId.size(); n++)
//     {
//         for(unsigned int m = 0; m < tracesId.size(); m++)
//         {
//             // Chiedo che n sia diverso da m, in questo modo avro tutte i punti due volte ma nelle funzioni successive mi facilita la ricerca
//             if(n != m)
//             {
//                 Vector3d vec = trace.coordinates[n].col(1) - trace.coordinates[n].col(0);
//                 Vector3d point = trace.coordinates[n].col(0);
//                 Vector3d vec1 = trace.coordinates[m].col(1) - trace.coordinates[m].col(0);
//                 Vector3d point1 = trace.coordinates[m].col(0);
//                 if((vec1.cross(vec)).norm() != 0)
//                 {
//                     double t;
//                     double u;
//                     t = ((point1.cross(vec1) - point.cross(vec1)).dot(vec.cross(vec1))) / (((vec.cross(vec1)).norm())*((vec.cross(vec1)).norm()));
//                     u = ((point.cross(vec) - point1.cross(vec)).dot(vec1.cross(vec))) / (((vec1.cross(vec)).norm())*((vec1.cross(vec)).norm()));

//                     Vector3d intersection = point + t * vec;
//                     Vector3d verify = point1 + u * vec1;

//                     // Se il punto di intersezione coincide con il punto di verifica, allora è corretto, lo inserisco in un vettore nella posizione corrispondente alla traccia intersecata dalla traccia che sto esaminando
//                     if((intersection - verify).norm() < tol && u >= 0 && u <= 1)
//                     {
//                         tuple<Vector3d, unsigned int, unsigned int> intersectionPoint;
//                         intersectionPoint = make_tuple(intersection, n, m);
//                         traceIntersection.push_back(intersectionPoint);
//                     }
//                 }
//             }
//         }
//     }
//     return traceIntersection;
// }


// bool notPassingCut(vector<vector<Vector3d>>& subfractureVertices, FractureNetwork::Traces trace, unsigned int traceId, vector<tuple<Vector3d, unsigned int, unsigned int>> tracesIntersection, vector<vector<Vector3d>>& subfractureVertices1)
// {
//     vector<Vector3d> traceCoord;
//     traceCoord.push_back(trace.coordinates[traceId].col(0));
//     traceCoord.push_back(trace.coordinates[traceId].col(1));
//     vector<tuple<Vector3d, unsigned int, unsigned int>> intersections;
//     // Salvo solo i punti di intersezione tra la traccia che sto studiando e tutte le altre
//     for(unsigned int i = 0; i < tracesIntersection.size(); i++)
//     {
//         if(get<1>(tracesIntersection[i]) == traceId)
//             intersections.push_back(tracesIntersection[i]);
//     }
//     // Giro su tutte le sottofratture che ho trovato fino a questo momento
//     for(unsigned int n = 0; n < subfractureVertices.size(); n++)
//     {
//         vector<Vector3d> subfracture = subfractureVertices[n];
//         unsigned int c = 0;
//         vector<Vector3d> subtrace;
//         // Giro su tutti i punti di intersezione tra la traccia che sto studiando e tutte le altre
//         for(unsigned int m = 0; m < intersections.size(); m++)
//         {
//             Vector3d point = get<0>(intersections[m]);
//             // Giro su tutti i lati della sottofrattura
//             for(unsigned int s = 0; s < subfracture.size(); s++)
//             {
//                 unsigned int t = (s + 1) % subfracture.size();
//                 Vector3d vec = subfracture[t] - subfracture[s];
//                 Vector3d P = subfracture[s];
//                 // Voglio vedere se il punto point si trova sui lati della sottofrattura
//                 double k = (vec.dot(point - P)) / (vec.norm() * vec.norm());
//                 if(k >= 0 && k <= 1)
//                 {
//                     c++;
//                     subtrace.push_back(point);
//                 }
//             }
//         }
//         // Se c = 2 vuol dire che ho trovato due punti di intersezione tra le tracce che appartengono alla sottofrattura che sto studiando, dunque posso consierare la
//         // sottotraccia come passante per questa sottofrattura e calcolo le due sottofratture
//         if(c == 2)
//         {
//             if(!createSubfracture(subfracture, subtrace, subfractureVertices1))
//             {
//                 cerr << "Error creating subtraces" << endl;
//                 return false;
//             }
//         }
//         // Se c = 1 vuol dire che trovo un solo punto di intersezione tra la mia traccia e le altre tracce che appartiene alla sottofrattura, dunque devo estendere la sottotraccia
//         else if(c == 1)
//         {
//             // DEVO METTERE DENTRO SUBTRACE UNO DEGLI ESTREMI DELLA TRACCIA, QUALE??? QUELLO CON DISTANZA MINIMA DALL'ALTRO PUNTO CHE HO IN SUBTRACE??
//             vector<Vector3d> extendedVertices = extendTraces(subfracture, subtrace);
//         }

//         // Se c = 0 vuol dire che non c'è nessun punto di intersezione tra la mia traccia e le altre che appartengono alla sottofrattura, dunque per questa
//     }


//     return true;
// }

// // Questa funzione verifica se un punto è interno ad un poligono
// bool pointInsidePolygon(vector<Vector3d> subFracture, Vector3d point)
// {
//     double tol = 1e-6;
//     unsigned int i = 0;
//     Vector3d vec = subFracture[i] - subFracture[i + 1];
//     Vector3d vec1 = subFracture[i + 2] - subFracture[i + 1];
//     Vector3d normal = vec.cross(vec1);
//     Vector3d v = subFracture[i + 1] - point;
//     double angleSum = 0;
//     // Verifico che il punto appartenga al piano della sottofrattura
//     if((v.cross(normal)).norm() < tol)
//     {
//         for(unsigned int n = 0; n < subFracture.size(); n++)
//         {
//             unsigned int k = (n+1) % subFracture.size();
//             Vector3d a = subFracture[n] - point;
//             Vector3d b = subFracture[k] - point;
//             double angle = acos(a.dot(b) / (a.norm() * b.norm()));
//             angleSum += angle;
//         }
//     }
//     else
//         return false;

//     // Verifico che l'angolo sia 2pi e quindi che il punto è interno al poligono
//     if(angleSum - (2 * M_PI) > tol)
//         return false;

//     return true;
// }


}

