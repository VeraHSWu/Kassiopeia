#include "KGAxialMesher.hh"
#include "KGCoreMessage.hh"
#include "KGInterfaceBuilder.hh"
#include "KGVTKAxialMeshPainter.hh"
#include "KXMLInitializer.hh"
#include "KXMLTokenizer.hh"

using namespace KGeoBag;
using namespace katrin;
using namespace std;

int main(int argc, char** argv)
{
    if (argc < 3) {
        cout << "usage: ./AxialMeshViewer <config_file_name.xml> <geometry_path> [...]" << endl;
        return -1;
    }

    coremsg(eNormal) << "starting initialization..." << eom;

    auto& tXML = KXMLInitializer::GetInstance();
    tXML.AddDefaultIncludePath(CONFIG_DEFAULT_DIR);
    tXML.Configure(argc, argv);

    deque<string> tPathList = tXML.GetArguments().ParameterList();
    tPathList.pop_front();  // strip off config file name

    coremsg(eNormal) << "...initialization finished" << eom;

    KVTKWindow tWindow;
    tWindow.SetName("KGeoBag Axial Mesh Viewer");
    tWindow.SetFrameColorRed(0.);
    tWindow.SetFrameColorGreen(0.);
    tWindow.SetFrameColorBlue(0.);
    tWindow.SetDisplayMode(true);
    tWindow.SetWriteMode(true);

    KGAxialMesher tMesher;

    KGVTKAxialMeshPainter tPainter;
    tPainter.SetName("AxialMeshPainter");
    tPainter.SetDisplayMode(true);
    tPainter.SetWriteMode(true);
    tPainter.SetColorMode(KGVTKAxialMeshPainter::sArea);

    for (auto& tPath : tPathList) {
        for (auto& tSurface : KGInterface::GetInstance()->RetrieveSurfaces(tPath)) {
            tSurface->MakeExtension<KGAxialMesh>();
            tSurface->AcceptNode(&tMesher);
            tSurface->AcceptNode(&tPainter);
        }
        for (auto& tSpace : KGInterface::GetInstance()->RetrieveSpaces(tPath)) {
            tSpace->MakeExtension<KGAxialMesh>();
            tSpace->AcceptNode(&tMesher);
            tSpace->AcceptNode(&tPainter);
        }
    }

    tWindow.AddPainter(&tPainter);
    tWindow.Render();
    tWindow.Write();
    tWindow.Display();
    tWindow.RemovePainter(&tPainter);

    return 0;
}
