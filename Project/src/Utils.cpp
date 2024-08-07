#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <string>
#include <limits>

using namespace std;

double tol = 1e-10;

// Questa funzione importa un file, studia le fratture al suo interno e stampa i risultati
bool ImportFracture(const string fileNameInput, const string fileNameOutput, const string fileNameOutputReordered,
                    const string filePathInput, const string filePathOutput , DiscreteFractureNetwork& fracture, Traces& trace)
{
    // Chiamo la funzione che legge il contenuto del file
    if(!ReadFracture(filePathInput, fileNameInput, fracture))
    {
        cerr << "Something wrong with the reading of the fracture" << endl;
        return false;
    }
    // Chiamo la funzione che trova le intersezioni tra le fratture
    if(!FractureIntersection(fracture, trace))
    {
        cerr << "Something wrong while computing fracture intersection" << endl;
        return false;
    }
    // Chiamo la funzione che stampa le tracce su file
    if(!PrintOnFile(fileNameOutput, filePathOutput, trace))
    {
        cerr << "Something wrong while printing the result in a file" << endl;
        return false;
    }
    // Chiamo la funzione che riordina le tracce per lunghezza
    if(!TraceReorder(fracture, trace))
    {
        cerr << "Something wrong while reordering traces" << endl;
        return false;
    }
    // Chiamo la funzione che stampa le tracce riordinate
    if(!printTraces(fileNameOutputReordered, filePathOutput, trace, fracture))
    {
        cerr << "Something wrong printing the traces reordered" << endl;
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


// Questa funzione crea una bounding box attorno alla frattura
vector<Vector3d> BBox3D(const MatrixXd& vertices)
{
    // Creo il vettore per la bounding box con due elementi, il minimo e il massimo
    vector<Vector3d> bbox(2);
    //Inizializzato a valori infinitamente grandi positivi e negativi
    bbox[0] = Vector3d::Constant(numeric_limits<double>::infinity());
    bbox[1] = Vector3d::Constant(-numeric_limits<double>::infinity());
    //Aggiorna quando trova un nuovo massimo o minimo elemento per elemento
    for (int i = 0; i < vertices.cols(); ++i) {
        bbox[0] = bbox[0].cwiseMin(vertices.col(i));
        bbox[1] = bbox[1].cwiseMax(vertices.col(i));
    }
    return bbox;
}


// Questa funzione calcola l'intersezione tra due piani creati da due fratture
bool FractureIntersection(const DiscreteFractureNetwork fracture, Traces& trace)
{
    // Creo il vettore per la bounding box
    vector<tuple<unsigned int, vector<Vector3d>>> BBoxVect;
    // Calcolo la bounding box per ogni frattura presente nel file
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        unsigned int ID = fracture.fractureID[i];
        vector<Vector3d> bbox = BBox3D(fracture.vertices[ID]);
        BBoxVect.push_back(make_tuple(ID, bbox));
    }
    // Presa una frattura
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        // Ciclo su tutte le fratture successive
        for(unsigned int j = 0; j < fracture.numFracture; j++)
        {
            if(i < j)
            {
                // Verifico se c'è intersezione tra le bounding box delle due fratture
                bool intersect = !(get<1>(BBoxVect[i])[0].x() > get<1>(BBoxVect[j])[1].x() ||
                                      get<1>(BBoxVect[i])[0].y() > get<1>(BBoxVect[j])[1].y() ||
                                      get<1>(BBoxVect[i])[0].y() > get<1>(BBoxVect[j])[1].z() ||
                                      get<1>(BBoxVect[j])[1].x() < get<1>(BBoxVect[i])[0].x() ||
                                      get<1>(BBoxVect[j])[1].y() < get<1>(BBoxVect[i])[0].y() ||
                                      get<1>(BBoxVect[j])[1].z() < get<1>(BBoxVect[i])[0].z());
                // Se non c'è intersezione salto tutto il codice successivo e prendo la frattura successiva
                if (!intersect)
                    continue;

                // Prendo gli id delle due fratture
                unsigned int Id1 = get<0>(BBoxVect[i]);
                unsigned int Id2 = get<0>(BBoxVect[j]);

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

                // Calcolo i termini noti
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


// Questa funzione trova le tracce come intersezioni di due fratture
bool FindTraces(const Vector3d s, const Vector3d point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2, Traces& trace)
{
    unsigned int n = 0;
    unsigned int m = 0;
    vector<Vector3d> Point1;
    vector<Vector3d> Point2;

    // Nel ciclo while prendo ogni volta un segmento della prima e un segmento della seconda frattura
    while(n < fracture.NumVertices[Id1] && m < fracture.NumVertices[Id2])
    {
        // Calcolo gli indici del vertice successivo al vertice considerato, se considero l'ultimo vertice il successivo sarà lo 0.
        unsigned int k = (n+1) % fracture.NumVertices[Id1];
        unsigned int h = (m+1) % fracture.NumVertices[Id2];

        // Prendo i vettori tra i due vertici
        Vector3d v1;
        v1 = fracture.vertices[Id1].col(n) - fracture.vertices[Id1].col(k);
        Vector3d P1;
        P1 = fracture.vertices[Id1].col(k);
        Vector3d v2;
        v2 = fracture.vertices[Id2].col(m) - fracture.vertices[Id2].col(h);
        Vector3d P2;
        P2 = fracture.vertices[Id2].col(h);

        // Verifico che le rette non siano parallele alla retta di intersezione tra i piani
        if((v1.cross(s)).norm() > tol * max(v1.norm(), s.norm()) || (v2.cross(s)).norm() > tol * max(v2.norm(), s.norm()))
        {
            double t1;
            double u1;
            t1 = ((P1.cross(v1) - point.cross(v1)).dot(s.cross(v1))) / (s.cross(v1).squaredNorm());
            u1 = ((point.cross(s) - P1.cross(s)).dot(v1.cross(s))) / (v1.cross(s).squaredNorm());
            double t2;
            double u2;
            t2 = ((P2.cross(v2) - point.cross(v2)).dot(s.cross(v2))) /  (s.cross(v2).squaredNorm());
            u2 = ((point.cross(s) - P2.cross(s)).dot(v2.cross(s))) / (v2.cross(s).squaredNorm());

            // Calcolo i punti di intersezione e i punti di verifica
            Vector3d intersection1 = point + t1 * s;
            Vector3d verify1 = P1 + u1 * v1;
            Vector3d intersection2 = point + t2 * s;
            Vector3d verify2 = P2 + u2 * v2;

            // Se il punto di intersezione coincide con il punto di verifica, allora è corretto
            if((intersection1 - verify1).norm() < tol * max(intersection1.norm(), verify1.norm()) && u1 >= 0 && u1 <= 1)
                Point1.push_back(intersection1);

            if((intersection2 - verify2).norm() < tol * max(intersection2.norm(), verify2.norm()) && u2 >= 0 && u2 <= 1)
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
        // Creo un vettore unico con i 4 punti di intersezione
        vector<Vector3d> intersectionPoints = {Point1[0], Point1[1], Point2[0], Point2[1]};

        // Ottengo le posizioni relative dei quattro punti sulla retta di intersezione tra i due piani, a e b saranno relativi alla prima frattura, c e d alla seconda
        double a = ((intersectionPoints[0] - point).dot(s)) / (s.norm() * s.norm());
        double b = ((intersectionPoints[1] - point).dot(s)) / (s.norm() * s.norm());
        double c = ((intersectionPoints[2] - point).dot(s)) / (s.norm() * s.norm());
        double d = ((intersectionPoints[3] - point).dot(s)) / (s.norm() * s.norm());


        // Confronto le posizioni relative dei punti per capire come sono posizionati e quali sono i due estremi della traccia, quindi salvo la traccia
        // Innanzitutto riordino le posizioni rispetto alla frattura, in modo da avere a prima di b e c prima di d
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
        // In questo caso c si trova tra a e b e b si trova tra c e d, dunque la traccia è quella formata da c e b
        if(c > a && c < b && b < d)
            SaveTraces(c, b, point, s, trace, Id1, Id2);
        // In questo caso a si trova tra c e d e b si trova dopo d, dunque la traccia è quella formata da a e d
        else if(a > c && a < d && d < b)
            SaveTraces(a, d, point, s, trace, Id1, Id2);
        // In questo caso c si trova tra a e d e b si trova dopo d, dunque la traccia è quella formata da c e d
        else if(c > a && c < d && d < b)
            SaveTraces(c, d, point, s, trace, Id1, Id2);
        // In questo caso a si trova tra c e b e d si trova dopo b, dunque la traccia è quella formata da a e b
        else if(a > c && a < b && b < d)
            SaveTraces(a, b, point, s, trace, Id1, Id2);
        // In questo caso a e c coincidono e b e d coincidono, dunque i due segmenti sono sovrapposti e la traccia è formata da una delle due coppie
        else if(fabs((a - c)) <= tol * max(fabs(a), fabs(c)) && fabs((b - d)) <= tol * max(fabs(b), fabs(d)))
            SaveTraces(a, b, point, s, trace, Id1, Id2);
        // In questo caso a e c coincidono mentre d si trova dopo b, dunque la traccia è formata da a e b
        else if(fabs((a - c)) <= tol * max(fabs(a), fabs(c)) && b < d)
            SaveTraces(a, b, point, s, trace, Id1, Id2);
        // In questo caso a e c coincidono mentre b si trova dopo d, dunque la traccia è formata da c e d
        else if(fabs((a - c)) <= tol * max(fabs(a), fabs(c)) && d < b)
            SaveTraces(c, d, point, s, trace, Id1, Id2);
        // In questo caso b e d coincidono mentre a si trova dopo c, dunque la traccia è formata da a e b
        else if(fabs((b - d)) <= tol * max(fabs(b), fabs(d)) && c < a)
            SaveTraces(a, b, point, s, trace, Id1, Id2);
        // In questo caso b e d coincidono menyte c si trova dopo a, dunque la traccia è formata da c e d
        else if(fabs((b - d)) <= tol * max(fabs(b), fabs(d)) && a < c)
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
    // Memorizzo gli estremi della traccia
    trace.coordinates.push_back(Coordinates);
    // Memorizzo gli id delle due fratture che generano la traccia
    Vector2i idFracture;
    idFracture << Id1, Id2;
    trace.fractureId.push_back(idFracture);
    // Memorizzo l'id della traccia
    trace.traceId.push_back(trace.numTraces);
    // Memorizzo la lunghezza della traccia
    double lengthSegment = sqrt(((P[0]-Q[0])*(P[0]-Q[0])+((P[1]-Q[1])*(P[1]-Q[1]))+((P[2]-Q[2])*(P[2]-Q[2]))));
    trace.length.push_back(lengthSegment);
    // Incremento il numero di tracce
    trace.numTraces++;
}


// Questa funzione stampa i risultati in un file
bool PrintOnFile(const string fileName, const string filePath, Traces trace)
{
    ofstream file;
    file.open(filePath + fileName);
    // Se fallisce l'apertura del file sollevo un errore
    if(file.fail())
    {
        cerr << "Error opening the file" << endl;
        return false;
    }
    // Prima stampo il numero di tracce
    file << "# Number of traces" << endl;
    file << trace.numTraces << endl;
    // Poi stampo le informazioni per ogni traccia, l'id, gli id delle due fratture che generano la traccia e le coordinate degli estremi della traccia
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


// Questa funzione riordina le tracce in passanti e non passanti e in ordine di lunghezza
bool TraceReorder(DiscreteFractureNetwork& fracture, Traces& trace)
{
    // Per ogni frattura
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        vector<tuple<unsigned int, bool, double>> fractureTraces;
        vector<unsigned int> passing = {};
        vector<unsigned int> notPassing = {};
        vector<double> lengthPassing = {};
        vector<double> lengthNotPassing = {};
        // Considero una traccia alla volta
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
                    if((PQ.cross(PA)).norm() < tol * max(PQ.norm(), PA.norm()) || (PQ.cross(PB)).norm() < tol * max(PQ.norm(), PB.norm()))
                        counter++;
                    n++;
                }
                // Finito il ciclo su tutti i segmenti, se il contatore è uguale a 2 vuol dire che i vertici della traccia appartengono a due segmenti della frattura,
                // dunque la traccia è passante per quella frattura. Memorizzo le posizioni, che corrispondono anche agli id delle tracce, in due vettori, uno per le tracce
                // passanti e uno per quelle non passanti.
                if(counter == 2)
                    passing.push_back(position);
                else
                    notPassing.push_back(position);
            }
        }
        // Creo altri due vettori, uno per la lunghezza delle tracce passanti e uno per la lunghezza delle tracce non passanti e inserisco le lunghezze mantenendo l'ordine degli id delle tracce.
        for(unsigned int a : passing)
            lengthPassing.push_back(trace.length[a]);
        for(unsigned int a : notPassing)
            lengthNotPassing.push_back(trace.length[a]);

        // Ora riordino il vettore delle lunghezze in ordine decrescente scambiando anche le posizioni del vettore di indici, così da mantenere la corrispondenza.
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
        // Creo le triplette (id, tips, lunghezza) come richiesto dall'output, prima per le passanti e poi per le non passanti e memorizzo le triplette
        // nel vettore fractureTraces
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
        // Aggiungo fractureTraces alla struttura
        trace.traceReordered.push_back(fractureTraces);
    }
    return true;
}


// Questa funzione riordina le tracce per lunghezza
bool reordering(vector<unsigned int>& idTraces, vector<double>& length)
{
    // Se ho un numero diverso di id e di lunghezze, allora restituisco false perchè c'è un errore
    if(idTraces.size() != length.size())
        return false;
    // Creo un vettore di coppie (lunghezza, id), una coppia per ogni traccia
    vector<pair<double, unsigned int>> couple;
    for(unsigned int i = 0; i < idTraces.size(); i++)
        couple.push_back((make_pair(length[i], idTraces[i])));
    // Riordino le coppie in ordine decrescente di lunghezza
    sort(couple.begin(), couple.end(), [](const pair<double, unsigned int>& a, const pair<double, unsigned int>& b){
        return a.first > b.first;
    });
    // Ora memorizzo l'ordine delle coppie nei due vettori iniziali degli id e delle lunghezze
    for(size_t i = 0; i < couple.size(); i++)
    {
        length[i] = couple[i].first;
        idTraces[i] = couple[i].second;
    }
    return true;
}


// Questa funzione stampa le tracce riordinate in un file
bool printTraces(const string fileName, const string filePath, Traces trace, DiscreteFractureNetwork fracture)
{
    ofstream file;
    file.open(filePath + fileName);
    // Se fallisce l'apertura del file di output viene sollevato un errore
    if(file.fail())
    {
        cerr << "Error opening the file" << endl;
        return false;
    }
    // Per ogni frattura
    for(unsigned int i = 0; i < fracture.numFracture; i++)
    {
        unsigned int fractureNumTraces = trace.traceReordered[i].size();
        // Salvo l'id della frattura e il numero di tracce presenti in questa frattura
        file << "# FractureId; NumTraces" << endl;
        file << fracture.fractureID[i] << "; " << fractureNumTraces << endl;
        // Poi salvo, per ogni traccia, il suo id, il booleano che mi indica se la traccia è passante o meno e la lunghezza della traccia
        file << "# TraceId; Tips; Length" << endl;
        if(fractureNumTraces != 0)
            for(unsigned int j = 0; j < fractureNumTraces; j++)
                file << get<0>(trace.traceReordered[i][j]) << "; " << get<1>(trace.traceReordered[i][j]) << "; " << get<2>(trace.traceReordered[i][j]) << endl;
        file << endl;
    }
    file.close();
    return true;
}


// **********************************************************************************************************************************************


// Questa funzione effettua il taglio delle fratture in sottopoligoni e salva i risultati in una mesh poligonale
bool fractureCut(DiscreteFractureNetwork& fracture, Traces& trace)
{
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

                // Ora posso dividere la frattura in due sottofratture
                if(!createSubfracture(fractureCoord, traceCoord, subfractureVertices1))
                {
                    cerr << "Error creating subtraces" << endl;
                    return false;
                }
                Vector3d cuttingTrace = traceCoord[1] - traceCoord[0];
                // Ora voglio dividere le sottofratture in destra e sinistra
                // Per ogni sottofrattura creata
                for(unsigned int h = 0; h < subfractureVertices1.size(); h++)
                {
                    Vector3d vec;
                    // Giro sui vertici della sottofrattura
                    for(unsigned int k = 0; k < subfractureVertices1[h].size(); k++)
                    {
                        // Entro nell'if solo se il punto che sto considerando, cioè quello con indice k, non è uno dei due estremi della traccia per cui ho tagliato.
                        if((subfractureVertices1[h][k] - traceCoord[0]).norm() > tol * max((subfractureVertices1[h][k]).norm(), (traceCoord[0]).norm())
                            && (subfractureVertices1[h][k] - traceCoord[1]).norm() > tol * max((subfractureVertices1[h][k]).norm(), (traceCoord[1]).norm()))
                        {
                            // Memorizzo il vettore ed esco dal ciclo
                            vec = subfractureVertices1[h][k] - traceCoord[0];
                            k = subfractureVertices1[h].size();
                        }
                    }

                    // Calcolo il prodotto vettoriale tra la traccia che ho tagliato e il vettore tra un estremo della traccia e un vertice della frattura
                    Vector3d crossProduct = cuttingTrace.cross(vec);
                    // Se il prodotto scalare tra la normale alla frattura e il prodotto vettoriale è positivo, vuol dire che la normale alla frattura
                    // e il prodotto vettoriale hanno lo stesso verso, dunque la sottofrattura sta a sinistra della traccia che ha tagliato
                    if(normal.dot(crossProduct) > 0)
                        subfractureVerticesSx = subfractureVertices1[h];
                    // Altrimenti i due vettori hanno verso opposto, dunque la sottofrattura sta a destra della traccia che ha tagliato
                    else
                        subfractureVerticesDx = subfractureVertices1[h];
                }
                // Una volta che ho diviso i suoi elementi in dx e sx, svuoto subfractureVertices1
                subfractureVertices1.clear();

                // Ora voglio dividere le tracce che devono ancora tagliare le sottofratture tra quelle che tagliano la sottofrattura a destra e
                // quelle che tagliano la sottofrattura a sinistra.
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
                        // Se i prodotti scalari tra la normale alla frattura e i due prodotti vettoriali sono entrambi positivi o nulli,
                        // significa che la traccia si trova a sinistra della traccia per cui ho tagliato
                        if(normal.dot(crossProduct1) >= 0 && normal.dot(crossProduct2) >= 0)
                            subTracesSx.push_back(fractureTraces[s]);
                        // Se i prodotti scalari tra la normale alla frattura e i due prodotti vettoriali sono entrambi negativi o nulli,
                        // significa che la traccia si trova a destra della traccia per cui ho tagliato
                        else if(normal.dot(crossProduct1) <= 0 && normal.dot(crossProduct2) <= 0)
                            subTracesDx.push_back(fractureTraces[s]);
                        // Se un prodotto scalare è positivo e l'altro è negativo vuol dire che le tracce si intersecano, trovo il punto
                        // di intersezione e divido le sottotracce in destra e sinistra
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
                            // Se il punto di intersezione coincide con il punto di verifica e t è compreso tra 0 e 1, allora è corretto e si trova sulla traccia, dunque lo salvo come punto di intersezione tra le due tracce
                            if((intersection - verify).norm() < tol * max(intersection.norm(), verify.norm()) && t >= 0 && t <= 1)
                                traceIntersection = intersection;
                            // La sottotraccia che sta a destra sarà quella che va dal punto di intersezione al punto il cui vettore aveva prodotto scalare negativo
                            vector<Vector3d> newTraceDx = {traceIntersection, get<0>(fractureTraces[s])[1]};
                            pair<vector<Vector3d>,unsigned int> traceDx = make_pair(newTraceDx, traceId1);
                            subTracesDx.push_back(traceDx);
                            // La sottotraccia che sta a sinistra sarà quella che va dal punto di intersezione al punto il cui vettore aveva prodotto scalare positivo
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
                            // Se il punto di intersezione coincide con il punto di verifica e t è compreso tra 0 e 1, allora è corretto e si trova sulla traccia, dunque lo salvo come punto di intersezione tra le due tracce
                            if((intersection - verify).norm() < tol * max(intersection.norm(), verify.norm()) && t >= 0 && t <= 1)
                                traceIntersection = intersection;
                            // La sottotraccia che sta a destra sarà quella che va dal punto di intersezione al punto il cui vettore aveva prodotto scalare negativo
                            vector<Vector3d> newTraceDx = {traceIntersection, get<0>(fractureTraces[s])[0]};
                            pair<vector<Vector3d>,unsigned int> traceDx = make_pair(newTraceDx, traceId1);
                            subTracesDx.push_back(traceDx);
                            // La sottotraccia che sta a sinistra sarà quella che va dal punto di intersezione al punto il cui vettore aveva prodotto scalare positivo
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

        // Ora salvo i risultati nella mesh
        PolygonalMesh mesh = {};
        if(!createMesh(subFracture, mesh))
        {
            cerr << "Error creating the mesh" << endl;
            return false;
        }
    }
    return true;
}


// Questa funzione effettua il taglio di una frattura per una traccia passante
bool createSubfracture(vector<Vector3d> subfracture, vector<Vector3d> cuttingTrace, vector<vector<Vector3d>>& subfractureVertices1)
{
    unsigned int numVertices = subfracture.size();
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
        Vector3d subFracVec = subfracture[k] - subfracture[n];
        // Creo un vettore dal punto della sottofrattura ad un estremo della traccia
        Vector3d firstVert = cuttingTrace[0] - subfracture[n];
        // E un vettore dal punto della frattura all'altro estremo della traccia
        Vector3d secondVert = cuttingTrace[1] - subfracture[n];
        // Inserisco il punto che sto considerando in Points
        Points.push_back(subfracture[n]);
        c++;
        // Se il prodotto vettoriale tra il vettore del lato della frattura e il vettore ad un estremo della traccia è 0, allora questo estremo
        // della traccia si trova sul lato della frattura, dunque lo salvo in Points e memorizzo la sua posizione all'interno di quest'ultimo
        if((subFracVec.cross(firstVert)).norm() < tol * max(subFracVec.norm(), firstVert.norm()))
        {
            Points.push_back(cuttingTrace[0]);
            pos1 = c;
            c++;
        }
        else if((subFracVec.cross(secondVert)).norm() < tol * max(subFracVec.norm(), secondVert.norm()))
        {
            Points.push_back(cuttingTrace[1]);
            pos2 = c;
            c++;
        }
        n++;
    }

    vector<Vector3d> subfracture1, subfracture2;
    unsigned int m = 0;
    // Con questo ciclo prendo gli elementi di Points e li salvo in subfracture1 finchè non trovo la posizione che ho memorizzato per un estremo della traccia.
    // Quando la trovo passo alla posizione dell'altro estremo della traccia e finisco il giro su Points. In questo modo salvo in subfracture1 gli elementi della
    // prima sottofrattura e "taglio fuori" quelli della seconda.
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
    // Ora devo salvare gli elementi della seconda sottofrattura, se l'estremo memorizzato in pos1 viene prima di quello in pos2, parto da pos1
    // e memorizzo i punti in subfracture2 finchè non trovo pos2.
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
    // Altimenti parto da pos2 e salvo gli elementi finchè non trovo pos1.
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
    // Memorizzo le due sottofratture che ho trovato in un unico vettore.
    subfractureVertices1.push_back(subfracture1);
    subfractureVertices1.push_back(subfracture2);

    return true;
}


// Questa funzione estende le tracce non passanti
vector<Vector3d> extendTraces(vector<Vector3d> subFractureVertices, vector<Vector3d> subTraceVertices)
{
    // vec è il vettore della sottotraccia da estendere
    Vector3d vec = subTraceVertices[1] - subTraceVertices[0];
    Vector3d point = subTraceVertices[0];
    vector<Vector3d> extendedVertices;
    unsigned int i = 0;
    // Per ogni vertice della sottofrattura
    while(i < subFractureVertices.size())
    {
        unsigned int k = (i + 1) % subFractureVertices.size();
        Vector3d vec1 = subFractureVertices[k] - subFractureVertices[i];
        Vector3d point1 = subFractureVertices[i];
        // Verifico che la traccia da estendere e il lato della sottofrattura che sto considerando non siano paralleli, se lo fossero passo direttamente
        // al lato successivo della sottofrattura
        if((vec1.cross(vec)).norm() > tol * max(vec1.norm(), vec.norm()))
        {
            double t;
            double u;
            t = ((point1.cross(vec1) - point.cross(vec1)).dot(vec.cross(vec1))) / (((vec.cross(vec1)).norm())*((vec.cross(vec1)).norm()));
            u = ((point.cross(vec) - point1.cross(vec)).dot(vec1.cross(vec))) / (((vec1.cross(vec)).norm())*((vec1.cross(vec)).norm()));

            // Trovo il punto di intersezione tra la retta che passa per il segmento e la retta che passa per il lato della sottofrattura che sto considerando
            Vector3d intersection = point + t * vec;
            Vector3d verify = point1 + u * vec1;

            // Se il punto di intersezione coincide con il punto di verifica e u è compreso tra 0 e 1, allora il punto è corretto e si trova sul lato della sottofrattura
            if((intersection - verify).norm() < tol * max(intersection.norm(), verify.norm()) && u >= 0 && u <= 1)
            {
                // Inserisco questo punto come estremo della traccia estesa
                extendedVertices.push_back(intersection);
            }
        }
        i++;
    }
    return extendedVertices;
}


// Questa funzione salva i risultati su una mesh
bool createMesh(vector<pair<vector<Vector3d>, vector<pair<vector<Vector3d>, unsigned int>>>> subFracture, PolygonalMesh& mesh)
{
    // Riservo lo spazio che mi può servire sovrastimandone la dimensione
    unsigned int dimension = subFracture.size();
    mesh.coordinates0D.reserve(3*dimension);
    mesh.cellId0D.reserve(3*dimension);
    mesh.cellId1D.reserve(4*dimension);
    mesh.verticesId1D.reserve(4*dimension);
    mesh.cellId2D.reserve(dimension);
    mesh.numVertices2D.reserve(dimension);
    mesh.numEdges2D.reserve(dimension);
    mesh.verticesId2D.reserve(dimension);
    mesh.edgesId2D.reserve(dimension);

    unsigned int id0D = 0;
    unsigned int id1D = 0;
    unsigned int id2D = 0;
    // Per ogni sottofrattura
    for(unsigned int k = 0; k < subFracture.size(); k++)
    {
        vector<Vector3d> polygon = get<0>(subFracture[k]);
        vector<unsigned int> verticesId;
        vector<unsigned int> edgesId;
        // Considero un vertice alla volta
        for(unsigned int h = 0; h < polygon.size(); h++)
        {
            // Cell0D
            // E considero il vertice sucessivo del poligono
            unsigned int t = (h+1) % polygon.size();
            // Verifico che questi vertici non siano già presenti nella mesh
            auto begin0D = mesh.coordinates0D.begin();
            auto end0D = mesh.coordinates0D.end();
            auto it1 = find(begin0D, end0D, polygon[h]);
            auto it2 = find(begin0D, end0D, polygon[t]);
            int id1;
            int id2;

            // Se il vertice h non è presente nella mesh lo salvo e memorizzo il suo id
            if(it1 == end0D)
            {
                mesh.coordinates0D.push_back(polygon[h]);
                mesh.cellId0D.push_back(id0D);
                id1 = id0D;
                id0D++;
            }
            // Altrimenti cerco e memorizzo il suo id
            else
                id1 = distance(begin0D, it1);

            // Se il vertice t non è presente nella mesh lo salvo e memorizzo il suo id
            if(it2 == end0D)
            {
                mesh.coordinates0D.push_back(polygon[t]);
                mesh.cellId0D.push_back(id0D);
                id2 = id0D;
                id0D++;
            }
            // Altrimenti cerco e memorizzo il suo id
            else
                id2 = distance(begin0D, it2);

            verticesId.push_back(id1);

            // Ora voglio verificare che questo punto non sia interno ad un segmento, infatti se lo fosse devo cancellare il segmento che lo contiene e sostituirlo con due segmenti,
            // ognuno formato da un estremo del segmento precedente e questo punto.
            // Un esempio di questa casistica è la frattura 0 del file FR3, dove ho un'intersezione delle due tracce "a T".
            for(unsigned int p = 0; p < mesh.verticesId1D.size(); p++)
            {
                int firstId = mesh.verticesId1D[p][0];
                int secondId = mesh.verticesId1D[p][1];
                int coupleId = mesh.cellId1D[p];
                Vector3d vec = mesh.coordinates0D[secondId] - mesh.coordinates0D[firstId];
                Vector3d vec1 = mesh.coordinates0D[id1] - mesh.coordinates0D[firstId];
                Vector3d vec2 = mesh.coordinates0D[id1] - mesh.coordinates0D[secondId];
                // Se il punto id1 si trova sul segmento, elimino il segmento e lo sostituisco con i due sotto-segmenti formati dal punto id1
                if((vec.cross(vec1)).norm() < tol * max(vec.norm(), vec1.norm()) && vec.dot(vec1) > 0 && (-vec).dot(vec2) > 0)
                {
                    auto it = remove(mesh.verticesId1D.begin(), mesh.verticesId1D.end(), mesh.verticesId1D[p]);
                    mesh.verticesId1D.erase(it, mesh.verticesId1D.end());
                    Vector2i newCouple1 = {firstId, id1};
                    Vector2i newCouple2 = {id1, secondId};
                    mesh.verticesId1D.insert(it, newCouple1);
                    mesh.verticesId1D.push_back(newCouple2);
                    mesh.cellId1D.push_back(id1D);

                    // Ora aggiorno anche la Cell2D
                    for(unsigned int q = 0; q < mesh.edgesId2D.size(); q++)
                    {
                        // Aggiungo l'id del nuovo lato al vettore di lati
                        auto it = find(mesh.edgesId2D[q].begin(), mesh.edgesId2D[q].end(), coupleId);
                        if(it != mesh.edgesId2D[q].end())
                        {
                            mesh.edgesId2D[q].insert(next(it), id1D);
                            mesh.numEdges2D[q]++;
                        }
                        // Aggiungo l'id del nuovo vertice al vettore di vertici
                        auto itFirst  = find(mesh.verticesId2D[q].begin(), mesh.verticesId2D[q].end(), firstId);
                        auto itSecond  = find(mesh.verticesId2D[q].begin(), mesh.verticesId2D[q].end(), secondId);
                        // Se il vettore r contiene entrambi i punti, aggiungo il punto in mezzo
                        if(itFirst != mesh.verticesId2D[q].end() && itSecond != mesh.verticesId2D[q].end())
                        {
                            mesh.verticesId2D[q].insert(next(itFirst), id1);
                            mesh.numVertices2D[q]++;
                        }
                    }
                    id1D++;
                }
            }

            // Cell1D
            Vector2i vertices1D = {id1, id2};
            auto begin1D = mesh.verticesId1D.begin();
            auto end1D = mesh.verticesId1D.end();
            auto itCouple1 = find(begin1D, end1D, vertices1D);
            int idEdge;

            // Verifico se questa coppia è già presente nella mesh
            if(itCouple1 == end1D)
            {
                Vector2i vertices1Dswap = {id2, id1};
                auto itCouple2 = find(begin1D, end1D, vertices1Dswap);
                // Verifico che anche la coppia opposta non è presente nella mesh
                if(itCouple2 == end1D)
                {
                    // Se entrambe non sono presenti, allora posso salvarla come nuova coppia
                    mesh.verticesId1D.push_back(vertices1D);
                    mesh.cellId1D.push_back(id1D);
                    idEdge = id1D;
                    id1D++;
                }
                // Se invece la coppia scambiata è presente nella mesh salvo l'id a quel lato
                else
                    idEdge = distance(begin1D, itCouple2);
            }
            // Se la coppia è presente nella mesh salvo l'id a quel lato
            else
                idEdge = distance(begin1D, itCouple1);

            edgesId.push_back(idEdge);
        }

        // Cell2D
        mesh.cellId2D.push_back(id2D);
        id2D++;
        mesh.verticesId2D.push_back(verticesId);
        mesh.numVertices2D.push_back(verticesId.size());
        mesh.edgesId2D.push_back(edgesId);
        mesh.numEdges2D.push_back(edgesId.size());
    }
    mesh.numCell0D = mesh.cellId0D.size();
    mesh.numCell1D = mesh.cellId1D.size();
    mesh.numCell2D = mesh.cellId2D.size();

    return true;
}







