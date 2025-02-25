#include "KMessage.h"
#include "KSFieldFinder.h"
#include "KSMainMessage.h"
#include "KSRootMagneticField.h"
#include "KTextFile.h"
#include "KThreeVector.hh"
#include "KXMLInitializer.hh"
#include "KXMLTokenizer.hh"


using namespace Kassiopeia;
using namespace katrin;
using namespace KGeoBag;
using namespace std;

int main(int argc, char** argv)
{
    if (argc < 6) {
        cout
            << "usage: ./SimpleMagneticFieldCalculatorOverXYplane <config_file.xml> <z> <rmax> <dr> <output_file.txt> <magnetic_field_name1> [<magnetic_field_name2> <...>] "
            << endl;
        // output_file can be "-" (-> write to terminal)
        exit(-1);
    }

    mainmsg(eNormal) << "starting initialization..." << eom;

    auto& tXML = KXMLInitializer::GetInstance();
    tXML.AddDefaultIncludePath(CONFIG_DEFAULT_DIR);
    tXML.Configure(argc, argv, true);

    deque<string> tParameters = tXML.GetArguments().ParameterList();
    tParameters.pop_front();  // strip off config file name

    double Z = stod(tParameters[0]);
    double Rmax = stod(tParameters[1]);
    double dR = stod(tParameters[2]);

    string outFileName(tParameters[3]);
    ofstream outFileStream;
    streambuf* outFileBuf;
    if (outFileName == "-") {
        outFileBuf = cout.rdbuf();
    }
    else {
        outFileStream.open(outFileName);
        outFileBuf = outFileStream.rdbuf();
    }
    ostream outFile(outFileBuf);

    mainmsg(eNormal) << "...initialization finished" << eom;

    // initialize magnetic field
    KSRootMagneticField tRootMagneticField;

    for (size_t tIndex = 4; tIndex < tParameters.size(); tIndex++) {
        KSMagneticField* tMagneticFieldObject = getMagneticField(tParameters[tIndex]);
        tMagneticFieldObject->Initialize();
        tRootMagneticField.AddMagneticField(tMagneticFieldObject);
    }

    KThreeVector tMagneticField;

    for (double x = -Rmax; x <= Rmax; x += dR) {
        for (double y = -Rmax; y <= Rmax; y += dR) {

            KThreeVector tPosition(x, y, Z);

            try {
                tRootMagneticField.CalculateField(tPosition, 0.0, tMagneticField);
            }
            catch (...) {
                mainmsg(eWarning) << "error - cannot calculate field at position <" << tPosition << ">" << eom;
                continue;
            }

            mainmsg(eNormal) << "Magnetic Field at position " << tPosition << " is " << tMagneticField << ", magnitude " << tMagneticField.Magnitude() << eom;

            outFile << std::fixed << std::setprecision(16) << tPosition.X() << "\t" << tPosition.Y() << "\t"
                    << tPosition.Z() << "\t" << tMagneticField.X() << "\t" << tMagneticField.Y() << "\t"
                    << tMagneticField.Z() << endl;
        }
    }

    for (size_t tIndex = 4; tIndex < tParameters.size(); tIndex++) {
        KSMagneticField* tMagneticFieldObject = getMagneticField(tParameters[tIndex]);
        tMagneticFieldObject->Deinitialize();
    }

    if (outFileStream.is_open())
        outFileStream.close();

    return 0;
}
