#include "UCDUtilities.hpp"

int main()
{
    Gedim::UCDUtilities ucdUtils;
    ucdUtils.ExportUCDAscii(points, pointProperties, cells, cellProperties, "output.icp");

    return 0;

}
