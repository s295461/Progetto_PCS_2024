#include <iostream>
#include "Eigen/Eigen"
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"
#include "PolygonalMesh.hpp"


using namespace std;
using namespace Eigen;



int main()
{
    // Creo le strutture per le fratture e le tracce
    DiscreteFractureNetwork fracture;
    Traces trace;

    string filePathInput = "DFN";
    string filePathOutput = "Result";

    /// FR3
    string fileNameFR3 = "/FR3_data.txt";
    string fileNameOutputFR3 = "/FR3_traces.txt";
    string fileNameOutputReorderedFR3 = "/FR3_traces_reordered.txt";

    // Chiamo la funzione che importa il file, salva le fratture, trova le tracce e le riordina per lunghezza
    if(!ImportFracture(fileNameFR3, fileNameOutputFR3, fileNameOutputReorderedFR3, filePathInput, filePathOutput, fracture, trace))
        return 1;
    // Chiamo la funzione che taglia le fratture e crea la mesh
    if(!fractureCut(fracture, trace))
        return 1;
    // Chiamo le funzioni che svuotano le strutture frature e trace prima di studiare il file successivo
    clearDiscreteFractureNetwork(fracture);
    clearTraces(trace);



    /// FR10
    string fileNameFR10 = "/FR10_data.txt";
    string fileNameOutputFR10 = "/FR10_traces.txt";
    string fileNameOutputReorderedFR10 = "/FR10_traces_reordered.txt";

    // Chiamo la funzione che importa il file, salva le fratture, trova le tracce e le riordina per lunghezza
    if(!ImportFracture(fileNameFR10, fileNameOutputFR10, fileNameOutputReorderedFR10, filePathInput, filePathOutput, fracture, trace))
        return 1;
    // Chiamo la funzione che taglia le fratture e crea la mesh
    if(!fractureCut(fracture, trace))
        return 1;
    // Chiamo le funzioni che svuotano le strutture frature e trace prima di studiare il file successivo
    clearDiscreteFractureNetwork(fracture);
    clearTraces(trace);



    /// FR50
    string fileNameFR50 = "/FR50_data.txt";
    string fileNameOutputFR50 = "/FR50_traces.txt";
    string fileNameOutputReorderedFR50 = "/FR50_traces_reordered.txt";

    // Chiamo la funzione che importa il file, salva le fratture, trova le tracce e le riordina per lunghezza
    if(!ImportFracture(fileNameFR50, fileNameOutputFR50, fileNameOutputReorderedFR50, filePathInput, filePathOutput, fracture, trace))
        return 1;
    // Chiamo la funzione che taglia le fratture e crea la mesh
    if(!fractureCut(fracture, trace))
        return 1;
    // Chiamo le funzioni che svuotano le strutture frature e trace prima di studiare il file successivo
    clearDiscreteFractureNetwork(fracture);
    clearTraces(trace);



    /// FR82
    string fileNameFR82 = "/FR82_data.txt";
    string fileNameOutputFR82 = "/FR82_traces.txt";
    string fileNameOutputReorderedFR82 = "/FR82_traces_reordered.txt";

    // Chiamo la funzione che importa il file, salva le fratture, trova le tracce e le riordina per lunghezza
    if(!ImportFracture(fileNameFR82, fileNameOutputFR82, fileNameOutputReorderedFR82, filePathInput, filePathOutput, fracture, trace))
        return 1;
    // Chiamo la funzione che taglia le fratture e crea la mesh
    if(!fractureCut(fracture, trace))
        return 1;
    // Chiamo le funzioni che svuotano le strutture frature e trace prima di studiare il file successivo
    clearDiscreteFractureNetwork(fracture);
    clearTraces(trace);



    /// FR200
    string fileNameFR200 = "/FR200_data.txt";
    string fileNameOutputFR200 = "/FR200_traces.txt";
    string fileNameOutputReorderedFR200 = "/FR200_traces_reordered.txt";

    // Chiamo la funzione che importa il file, salva le fratture, trova le tracce e le riordina per lunghezza
    if(!ImportFracture(fileNameFR200, fileNameOutputFR200, fileNameOutputReorderedFR200, filePathInput, filePathOutput, fracture, trace))
        return 1;
    // Chiamo la funzione che taglia le fratture e crea la mesh
    if(!fractureCut(fracture, trace))
        return 1;
    // Chiamo le funzioni che svuotano le strutture frature e trace prima di studiare il file successivo
    clearDiscreteFractureNetwork(fracture);
    clearTraces(trace);



    /// FR362
    string fileNameFR362 = "/FR362_data.txt";
    string fileNameOutputFR362 = "/FR362_traces.txt";
    string fileNameOutputReorderedFR362 = "/FR362_traces_reordered.txt";

    // Chiamo la funzione che importa il file, salva le fratture, trova le tracce e le riordina per lunghezza
    if(!ImportFracture(fileNameFR362, fileNameOutputFR362, fileNameOutputReorderedFR362, filePathInput, filePathOutput, fracture, trace))
        return 1;
    // Chiamo la funzione che taglia le fratture e crea la mesh
    if(!fractureCut(fracture, trace))
        return 1;
    // Chiamo le funzioni che svuotano le strutture frature e trace prima di studiare il file successivo
    clearDiscreteFractureNetwork(fracture);
    clearTraces(trace);

    return 0;



}
