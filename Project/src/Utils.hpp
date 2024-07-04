#pragma once

#include "DiscreteFractureNetwork.hpp"
#include "PolygonalMesh.hpp"

using namespace std;


// Questa funzione importa un file, studia le fratture al suo interno e stampa i risultati
bool ImportFracture(const string fileNameInput, const string fileNameOutput, const string fileNameOutputReordered,
                    const string filePathInput, const string filePathOutput, DiscreteFractureNetwork& fracture, Traces& trace);

// Questa funzione svuota tutti gli elementi all'interno della struttura DiscreteFractureNetwork
void clearDiscreteFractureNetwork(DiscreteFractureNetwork& fracture);

// Questa funzione svuota tutti gli elementi all'interno della struttura Traces
void clearTraces(Traces& trace);

// Questa funzione apre un file, legge il contenuto e lo salva in strutture dati adeguate
bool ReadFracture(const string& filePath, const string& fileName, DiscreteFractureNetwork& fracture);

// Questa funzione crea una bounding box attorno alla frattura
vector<Vector3d> BBox3D(const MatrixXd& vertices);

// Questa funzione calcola l'intersezione tra due piani creati da due fratture
bool FractureIntersection(const DiscreteFractureNetwork fracture, Traces& trace);

// Questa funzione trova le tracce come intersezioni di due fratture
bool FindTraces(const Vector3d s, const Vector3d Point, const DiscreteFractureNetwork& fracture, unsigned int Id1, unsigned int Id2, Traces& trace);

// Con questa funzione salvo tutti i valori legati ad una traccia nella struttura Traces.
void SaveTraces(double n, double m, Vector3d point, Vector3d s, Traces& trace, unsigned int Id1, unsigned int Id2);

// Questa funzione stampa i risultati in un file
bool PrintOnFile(const string fileName, const string filePath, Traces trace);

// Questa funzione riordina le tracce in passanti e non passanti e in ordine di lunghezza
bool TraceReorder(DiscreteFractureNetwork& fracture, Traces& trace);

// Questa funzione riordina le tracce per lunghezza
bool reordering(vector<unsigned int>& idTraces, vector<double>& length);

// Questa funzione stampa le tracce riordinate in un file
bool printTraces(const string fileName, const string filePath, Traces trace, DiscreteFractureNetwork fracture);


// *****************************************************************************************************************************************************


// Questa funzione effettua il taglio delle fratture in sottopoligoni e salva i risultati in una mesh poligonale.
bool fractureCut(DiscreteFractureNetwork& fracture, Traces& trace);

// Questa funzione effettua il taglio di una frattura per una traccia passante
bool createSubfracture(vector<Vector3d> subfracture, vector<Vector3d> cuttingTrace, vector<vector<Vector3d>>& subfractureVertices);

// Questa funzione estende le tracce non passanti
vector<Vector3d> extendTraces(vector<Vector3d> subFractureVertices, vector<Vector3d> subTraceVertices);

// Questa funzione salva i risultati su una mesh
bool createMesh(vector<pair<vector<Vector3d>, vector<pair<vector<Vector3d>, unsigned int>>>> subFracture, PolygonalMesh& mesh);
