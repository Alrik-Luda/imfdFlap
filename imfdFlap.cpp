//imfdFlap, submitted as part of the assessment in the NTFD3 course at TU Bergakademie Freiberg by Alrik Luda 
//this file is meant to be compiled locally with g++
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <filesystem>
#include <cmath>
#include <sstream>

using namespace std;

namespace fs = std::filesystem;

void printContent(ofstream& file, vector<string>& content){
    for(const string &line : content){
        file << line << endl;
    }

    file << "// ************************************************************************* //";
}

void printHeader(ofstream& file, string& filename, string& location, string& className){

    string objectName = "    object      " + filename;
    className = "    class       " + className + ";";
    location = "    location    \"" + location + "\";";

    objectName = objectName + ";";

    vector<string> header={
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        "| =========                 |                                                 |",
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |",
        "|  \\\\    /   O peration     | Version:  2312 used by Alrik to make this       |",
        "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |",
        "|    \\\\/     M anipulation  |                                                 |",
        "\\*---------------------------------------------------------------------------*/",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    arch        \"LSB;label=32;scalar=64\";",
        location,
        className,
        objectName,
        "}",
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
        };

    for(const string &line : header){
        file << line << endl;
    }
}

void writeU(const string workDir, const vector<double>U,const vector<double>geometryProperties, const vector<string>dirNames){
    for(int i = 0; i<U.size();i++){
        
        //write double to stringstream and then to string
        ostringstream strsU;
        strsU.clear();
        strsU << "#calc \"" << U.at(i) << "*" << 1e-3*geometryProperties.at(6) << "*" << 1e-3*geometryProperties.at(2); 
        string valueU;
        valueU = strsU.str();

        string filename = "U";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "volVectorField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 1 -1 0 0 0 0];",
            " ",
            "internalField\tuniform (0 0 0);",
            " ",
            "boundaryField",
            "{",
                "\tinlet",
                "\t{",
                    "\t\ttype\tflowRateInletVelocity;",
                    "\t\tvolumetricFlowRate\ttable",
                        "\t\t\t(",
                            "\t\t\t\t(0 0)",
                            "\t\t\t\t(1 " + valueU + "\") // u*h*b",
                        "\t\t\t);",
                    "\t\tvalue\tuniform (0 0 0);",
                "\t}",
                "\toutlet",
                "\t{",
                    "\t\ttype\tzeroGradient;",
                "\t}",
                "\tflap",
                "\t{",
                    "\t\ttype\tmovingWallVelocity;",
                    "\t\tvalue\tuniform (0 0 0);",
                "\t}",
                "\tatmosphere",
                "\t{",
                    "\t\ttype\tslip;",
                "\t}",
                "\t\"(lowerWall|frontAndBack)\"",
                "\t{",
                    "\t\ttype\tnoSlip;",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeK(const string workDir, const vector<double>k,const vector<double>I, const vector<string>dirNames){
    for(int i = 0; i<k.size();i++){
        
        //write double to stringstream and then to string
        ostringstream strsK, strsI;
        strsK.clear();
        strsI.clear();
        strsK << "uniform " << k.at(i);
        strsI << I.at(i);
        string valueK;
        string valueI;
        valueK = strsK.str();
        valueI = strsI.str();

        string filename = "k";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "volScalarField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written using \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 2 -2 0 0 0 0];",
            " ",
            "internalField\t" + valueK + ";",
            " ",
            "boundaryField",
            "{",
                "\tinlet",
                "\t{",
                    "\t\ttype\tturbulentIntensityKineticEnergyInlet;",
                    "\t\tintensity\t" + valueI + ";",
                    "\t\tvalue\t" + valueK + ";",
                "\t}",
                "\toutlet",
                "\t{",
                    "\t\ttype\tfixedValue;",
                    "\t\tvalue\t$internalField;",
                "\t}",
                "\t\"(flap|atmosphere|lowerWall|frontAndBack)\"",
                "\t{",
                    "\t\ttype\tkqRWallFunction;",
                    "\t\tvalue\tuniform 1e-06;",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeOmega(const string workDir, const vector<double>omega, const vector<string>dirNames){
    for(int i = 0; i<omega.size();i++){
        
        //write double to stringstream and then to string
        ostringstream strsOmega;
        strsOmega.clear();
        strsOmega << "uniform " << omega.at(i);
        string valueOmega;
        valueOmega = strsOmega.str();

        string filename = "omega";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "volScalarField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 0 -1 0 0 0 0];",
            " ",
            "internalField\t" + valueOmega + ";",
            " ",
            "boundaryField",
            "{",
                "\tinlet",
                "\t{",
                    "\t\ttype\tturbulentMixingLengthFrequencyInlet;",
                    "\t\tmixingLength\t0.002; //5% channel height",
                    "\t\tk\tk;",
                    "\t\tvalue\t" + valueOmega + ";",
                "\t}",
                "\toutlet",
                "\t{",
                    "\t\ttype\tzeroGradient;",
                "\t}",
                "\t\"(flap|atmosphere|lowerWall|frontAndBack)\"",
                "\t{",
                    "\t\ttype\tomegaWallFunction;",
                    "\t\tvalue\t" + valueOmega + ";",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeP(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "p";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "volScalarField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 2 -2 0 0 0 0];",
            " ",
            "internalField\tuniform 0;",
            " ",
            "boundaryField",
            "{",
                "\tinlet",
                "\t{",
                    "\t\ttype\tzeroGradient;",
                "\t}",
                "\toutlet",
                "\t{",
                    "\t\ttype\tfixedMean;",
                    "\t\tmeanValue\tconstant 0;",
                    "\t\tvalue\tuniform 0;",
                "\t}",
                "\t\"(flap|atmosphere|lowerWall|frontAndBack)\"",
                "\t{",
                    "\t\ttype\tzeroGradient;",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeNut(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "nut";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "volScalarField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 2 -1 0 0 0 0];",
            " ",
            "internalField\tuniform 0;",
            " ",
            "boundaryField",
            "{",
                "\tinlet",
                "\t{",
                    "\t\ttype\tcalculated;",
                    "\t\tvalue\tuniform 0;",
                "\t}",
                "\toutlet",
                "\t{",
                    "\t\ttype\tcalculated;",
                    "\t\tvalue\tuniform 0;",
                "\t}",
                "\t\"(flap|atmosphere|lowerWall|frontAndBack)\"",
                "\t{",
                    "\t\ttype\tnutkWallFunction;",
                    "\t\tU\tU;",
                    "\t\tvalue\tuniform 0;",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writePhi(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "phi";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "surfaceScalarField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 3 -1 0 0 0 0];",
            " ",
            "internalField\tuniform 0;",
            " ",
            "boundaryField",
            "{",
                "\t\"(inlet|outlet)\"",
                "\t{",
                    "\t\ttype\tcalculated;",
                    "\t\tvalue\t$internalField;",
                "\t}",
                "\t\"(flap|atmosphere|lowerWall|frontAndBack)\"",
                "\t{",
                    "\t\ttype\tcalculated;",
                    "\t\tvalue\tuniform 0;",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writePointDisplacement(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "pointDisplacement";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "pointVectorField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 1 0 0 0 0 0];",
            " ",
            "internalField\tuniform (0 0 0);",
            " ",
            "boundaryField",
            "{",
                "\t\"(inlet|outlet)\"",
                "\t{",
                    "\t\ttype\tfixedValue;",
                    "\t\tvalue\tuniform (0 0 0);",
                "\t}",
                "\tflap",
                "\t{",
                    "\t\ttype\tfixedValue;",
                    "\t\tvalue\t$internalField;",
                "\t}",
                "\t\"(atmosphere|lowerWall|frontAndBack)\"",
                "\t{",
                    "\t\ttype\tslip;",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeD(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "D";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "volVectorField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[0 1 0 0 0 0 0];",
            " ",
            "internalField\tuniform (0 0 0);",
            " ",
            "boundaryField",
            "{",
                "\tflap",
                "\t{",
                    "\t\ttype\tsolidDisplacementFoamForce;",
                    "\t\tforceField\tsolidForce;",
                    "\t\tvalue\tuniform (0 0 0);",
                "\t}",
                "\tbottom",
                "\t{",
                    "\t\ttype\tfixedValue;",
                    "\t\tvalue\tuniform (0 0 0);",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeSolidForce(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "solidForce";
        string location = "0";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "volVectorField";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dimensions\t[1 2 2 0 0 0 0];",
            " ",
            "internalField\tuniform (0 0 0);",
            " ",
            "boundaryField",
            "{",
                "\tflap",
                "\t{",
                    "\t\ttype\tcalculated;",
                    "\t\tvalue\tuniform (0 0 0);",
                "\t}",
                "\tbottom",
                "\t{",
                    "\t\ttype\tcalculated;",
                    "\t\tvalue\tuniform (0 0 0);",
                "\t}",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeDynamicMeshDict(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "dynamicMeshDict";
        string location = "constant";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "dynamicFvMesh\tdynamicMotionSolverFvMesh;",
            " ",
            "motionSolverLibs\t(\"libfvMotionSolvers.so\");",
            " ",
            "motionSolver\tdisplacementLaplacian;",
            " ",
            "displacementLaplacianCoeffs {",
            "\tdiffusivity\tquadratic inverseDistance (flap); //more patches can be specified by changing the format to ... quadratic inverseDistance #(\"flap\",...,...) and adding the # of patches affected ",
            "}"
        };

        printContent(outfile,content);
    }
}

void writeTransportProperties(const string workDir, const vector<double>U, const vector<double>fluidProperties, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        ostringstream strsTP;
        strsTP.clear();
        strsTP << "nu\tnu [ 0 2 -1 0 0 0 0 ] " << fluidProperties.at(1)/fluidProperties.at(0) << ";";
        string valueTP;
        valueTP = strsTP.str();

        string filename = "transportProperties";
        string location = "constant";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "transportModel\tNewtonian;",
            " ",
            valueTP
        };

        printContent(outfile,content);
    }
}

void writeTurbulenceProperties(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "turbulenceProperties";
        string location = "constant";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //adjusting for laminar and turbulent flow
        if(i<3){
            vector<string> content = {
                "simulationType\tlaminar;",
            };
            printContent(outfile,content);
        }
        else{
            vector<string> content = {
                "simulationType\tRAS;",
                " ",
                "RAS",
                "{",
                    "\tRASModel\tkOmegaSST;",
                    " ",
                    "\tturbulence\ton;",
                    " ",
                    "\tprintCoeffs\ton;",
                "}"
            };
            printContent(outfile,content);
        };
    }
}

void writeThermalProperties(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        string filename = "thermalProperties";
        string location = "constant";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            "thermalStress\tno;"
        };

        printContent(outfile,content);
    }
}

void writeMechanicalProperties(const string workDir, const vector<double>U, const vector<double>solidProperties, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        ostringstream strsMP;
        strsMP.clear();
        strsMP << "rho" << endl << "{" << endl << "\ttype\tuniform;" << endl << "\tvalue\t" << solidProperties.at(0) << ";" << endl << "}" << endl << endl << "nu" << endl << "{" << endl << "\ttype\tuniform;" << endl << "\tvalue\t" << solidProperties.at(2) << ";" << endl << "}" << endl << endl << "E" << endl << "{" << "\ttype\tuniform;" << endl << "\tvalue\t" << solidProperties.at(1) << ";" << endl << "}";
        string valueMP;
        valueMP = strsMP.str();

        string filename = "mechanicalProperties";
        string location = "constant";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        //of course this could be compactly written with \n but this style is intended to allow maximum clarity and debugging
        vector<string> content = {
            valueMP,
            " ",
            "planeStress\tno; //yes in 2D cases or if the geometry is thin and the stress components acting out of the 2D plane are negligible",
            " ",
            "//planeStrain condition 2D strains are negligible, applicable to thick 3d solids"
        };

        printContent(outfile,content);
    }
}

void writeFluidPreciceDict(const string workDir, const vector<double>U, const vector<double>fluidProperties, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        ostringstream strsPD;
        strsPD.clear();
        strsPD << "\trho\trho [1 -3 0 0 0 0 0] " << fluidProperties.at(0) << ";";
        string valuePD;
        valuePD = strsPD.str();

        string filename = "preciceDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/fluid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        vector<string> content = {
                "preciceConfig\t\"../precice-config.xml\";",
                " ",
                "participant Fluid;",
                " ",
                "modules (FSI);",
                " ",
                "interfaces",
                "{",
                    "\tInterface1",
                    "\t{",
                        "\t\tmesh\tFluid-Mesh;",
                        "\t\tpatches\t(flap);",
                        "\t\tlocations\tfaceCenters;",
                        " ",
                        "\t\treadData",
                        "\t\t(",
                        "\t\t\tDisplacement",
                        "\t\t);",
                        " ",
                        "\t\twriteData",
                        "\t\t(",
                        "\t\t\tForce",
                        "\t\t);",
                    "\t};",
                "};",
                " ",
                "FSI",
                "{",
                valuePD + " //required for force calculation",
                "};"
            };
        
        printContent(outfile,content);
    }
}

void writeFluidControlDict(const string workDir, const vector<double>U, const vector<string>timestep, const vector<string>dirNames, string writeCompression){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        ostringstream strsCD;
        strsCD.clear();
        strsCD << "deltaT\t" << timestep.at(i) << ";";
        string valueCD;
        valueCD = strsCD.str();

        string filename = "controlDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/fluid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "application\tpimpleFoam;",
            " ",
            "startFrom\tstartTime;",
            " ",
            "startTime\t0;",
            " ",
            "stopAt\tendTime;",
            " ",
            "endTime\t10;",
            " ",
            valueCD,
            " ",
            "writeControl\tadjustableRunTime;",
            " ",
            "writeInterval\t1e-02;",
            " ",
            "purgeWrite\t0;",
            " ",
            "writeFormat\tascii;",
            " ",
            "writePrecision\t12;",
            " ",
            "writeCompression\t" + writeCompression + ";",
            " ",
            "timeFormat\tgeneral;",
            " ",
            "timePrecision\t6;",
            " ",
            "functions",
            "{",
                "\tpreCICE_Adapter",
                "\t{",
                    "\t\ttype\tpreciceAdapterFunctionObject;",
                    "\t\tlibs\t(\"libpreciceAdapterFunctionObject.so\");",
                "\t}",
            "}"
        };
        
        printContent(outfile,content);
    }
}

void writeFluidFvSolution(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        string filename = "fvSolution";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/fluid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "solvers",
            "{",
                "\t\"(p|pcorr|Phi)\"",
                "\t{",
                    "\t\tsolver\tGAMG;",
                    "\t\ttolerance\t1e-08;",
                    "\t\trelTol\t0.001;",
                    "\t\tminIter\t3;",
                    "\t\tsmoother\tDIC;",
                    "\t\tnPreSweeps\t1;",
                    "\t\tnPostSweeps\t2;",
                    "\t\tnFinestSweeps\t2;",
                    "\t\tscaleCorrection\ttrue;",
                    "\t\tdirectSolveCoarsestLevel\tfalse;",
                    "\t\tcacheAgglomeration\ton; //don't use with dynamic mesh refinement",
                    "\t\tnCellsInCoarsestLevel 850; //sqrt(no. of cells)",
                    "\t\tagglomerator\tfaceAreaPair;",
                    "\t\tmergeLevels\t1;",
                "\t}",
                " ",
                "\t\"(pFinal|pcorrFinal)\"",
                "\t{",
                    "\t\t$p;",
                    "\t\trelTol\t0;",
                "\t}",
                " ",
                "\tU",
                "\t{",
                    "\t\tsolver\tsmoothSolver;",
                    "\t\tsmoother\tsymGaussSeidel;",
                    "\t\ttolerance\t1e-6;",
                    "\t\trelTol\t1e-4;",
                    "\t\tminIter\t2;",
                "\t}",
                " ",
                "\tUFinal",
                "\t{",
                    "\t\t$U;",
                    "\t\trelTol\t0;",
                "\t}",
                " ",
                "\tcellDisplacement",
                "\t{",
                    "\t\tsolver\tsmoothSolver;",
                    "\t\tsmoother\tsymGaussSeidel;",
                    "\t\ttolerance\t1e-12;",
                    "\t\trelTol\t1e-2;",
                    "\t\tminIter\t2;",
                "\t}",
                " ",
                "\tcellDisplacementFinal",
                "\t{",
                    "\t\t$cellDisplacement;",
                    "\t\trelTol\t0;",
                "\t}",
                " ",
                "\t\"(k|omega|R|nuTilda)\"",
                "\t{",
                    "\t\tsolver\tPBiCGStab;",
                    "\t\tpreconditioner\tDILU;",
                    "\t\ttolerance\t1e-12;",
                    "\t\trelTol\t1e-2;",
                "\t}",
                " ",
                "\t\"(k|omega|R|nuTilda)Final\"",
                "\t{",
                    "\t\t\"$(k|omega|R|nuTilda)\";",
                    "\t\trelTol\t0;",
                "\t}",
            "}",
            " ",
            "PIMPLE",
            "{",
                "\tnCorrectors\t4;",
                " ",
                "\tnNonOrthogonalCorrectors\t1;",
                " ",
                "\tmomentumPredictor\ttrue;",
                " ",
                "\tnOuterCorrectors\t1;",
                " ",
                "\tconsistent\ttrue;",
                " ",
                "\tcorrectPhi\ttrue;",
                " ",
                "\trelaxationFactors",
                "\t{",
                    "\t\tequations",
                    "\t\t{",
                        "\t\t\t\".*\"\t1;",
                    "\t\t}",
                "\t}"
                " ",
                "\tresidualControl",
                "\t{",
                    "\t\tp",
                    "\t\t{",
                        "\t\t\trelTol\t0;",
                        "\t\t\ttolerance\t1e-3;",
                    "\t\t}"
                    " ",
                    "\t\t\"(U|k|p|omega)\"",
                    "\t\t{",
                        "\t\t\trelTol\t0;",
                        "\t\t\ttolerance\t1e-4;",
                    "\t\t}",
                "\t}",
            "}"
            " ",
            "potentialFlow",
            "{",
                "\tnNonOrthogonalCorrectors\t10;",
            "}"
        };
        
        printContent(outfile,content);
    }
}

void writeFluidFvSchemes(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        string filename = "fvSchemes";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/fluid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "ddtSchemes",
            "{",
                "\tdefault\tbackward;",
            "}",
            " ",
            "gradSchemes",
            "{",
                "\tdefault\tGauss linear;",
            "}",
            " ",
            "divSchemes",
            "{",
                "\tdefault\tnone;",
                "\tdiv(phi,U)\tGauss linearUpwind grad(U);",
                "\tdiv((nuEff*dev2(T(grad(U)))))\tGauss linear;",
                "\tdiv(phi,k)\tGauss linearUpwind grad(k);",
                "\tdiv(phi,omega)\tGauss linearUpwind grad(omega);",
                "\tdiv(phi,nuTilda)\tGauss linear;",
                "\tdiv(phi,R)\tGauss linear;",
                "\tdiv(R)\tGauss linear;",
            "}",
            " ",
            "laplacianSchemes",
            "{",
                "\tdefault\tGauss linear corrected;",
            "}",
            " ",
            "interpolationSchemes",
            "{",
                "\tdefault\tlinear;",
            "}",
            " ",
            "snGradSchemes",
            "{",
                "\tdefault\tcorrected;",
            "}",
            " ",
            "wallDist",
            "{",
                "\tmethod\tmeshWave;",
            "}",
            " ",
            "fluxRequired",
            "{",
                "\tdefault\tno;",
                "\tp;",
            "}"
        };
        
        printContent(outfile,content);
    }
}

void writeFluidDecomposeParDict(const string workDir, const vector<double>U, const vector<int>fluidCores, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsFDPD;
        strsFDPD.clear();
        strsFDPD << fluidCores.at(i);
        string valueFDPD;
        valueFDPD = strsFDPD.str();

        string filename = "decomposeParDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/fluid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "numberOfSubdomains\t" + valueFDPD + ";",
            " ",
            "method\tscotch;"
        };
        
        printContent(outfile,content);
    }
}

void writeFluidBlockMeshDict(const string workDir, const vector<double>U, const vector<double>geometryProperties, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        double channelLengthFraction = 0.2;

        ostringstream strsBMD;
        strsBMD.clear();

        strsBMD << "x0\t" << -1*1e-3*geometryProperties.at(0)*channelLengthFraction*0.5 << ";" << endl << "x1\t" << -1*1e-3*geometryProperties.at(3)*0.5 << ";" << endl << "x2\t" << 1e-3*geometryProperties.at(3)*0.5 << ";" << endl << "x3\t" << 1e-3*geometryProperties.at(0)*channelLengthFraction*0.5 << ";" << endl << endl << "y0\t0" << ";" << endl << "y1\t" << 1e-3*geometryProperties.at(4) << ";" << endl << "y2\t" << 1e-3*geometryProperties.at(6) << ";" << endl << endl << "z0\t" << -1*1e-3*geometryProperties.at(2)*0.5 << ";" << endl << "z1\t" << -1*1e-3*geometryProperties.at(5)*0.5 << ";" << endl << "z2\t" << 1e-3*geometryProperties.at(5)*0.5 << ";" << endl << "z3\t" << 1e-3*geometryProperties.at(2)*0.5 << ";";

        string valueBMD;
        valueBMD = strsBMD.str();


        string filename = "blockMeshDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/fluid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "//flap base centred at 0,0,0",
            " ",
            valueBMD,
            " ",
            "vertices",
            "(",
                "\t($x0 $y0 $z0)\t// 0",
                "\t($x1 $y0 $z0)\t// 1",
                "\t($x2 $y0 $z0)\t// 2",
                "\t($x3 $y0 $z0)\t// 3",
                " ",
                "\t($x0 $y1 $z0)\t// 4",
                "\t($x1 $y1 $z0)\t// 5",
                "\t($x2 $y1 $z0)\t// 6",
                "\t($x3 $y1 $z0)\t// 7",
                " ",
                "\t($x0 $y2 $z0)\t// 8",
                "\t($x1 $y2 $z0)\t// 9",
                "\t($x2 $y2 $z0)\t// 10",
                "\t($x3 $y2 $z0)\t// 11",
                " ",
                "\t($x0 $y0 $z1)\t// 12",
                "\t($x1 $y0 $z1)\t// 13",
                "\t($x2 $y0 $z1)\t// 14",
                "\t($x3 $y0 $z1)\t// 15",
                " ",
                "\t($x0 $y1 $z1)\t// 16",
                "\t($x1 $y1 $z1)\t// 17",
                "\t($x2 $y1 $z1)\t// 18",
                "\t($x3 $y1 $z1)\t// 19",
                " ",
                "\t($x0 $y2 $z1)\t// 20",
                "\t($x1 $y2 $z1)\t// 21",
                "\t($x2 $y2 $z1)\t// 22",
                "\t($x3 $y2 $z1)\t// 23",
                " ",
                "\t($x0 $y0 $z2)\t// 24",
                "\t($x1 $y0 $z2)\t// 25",
                "\t($x2 $y0 $z2)\t// 26",
                "\t($x3 $y0 $z2)\t// 27",
                " ",
                "\t($x0 $y1 $z2)\t// 28",
                "\t($x1 $y1 $z2)\t// 29",
                "\t($x2 $y1 $z2)\t// 30",
                "\t($x3 $y1 $z2)\t// 31",
                " ",
                "\t($x0 $y2 $z2)\t// 32",
                "\t($x1 $y2 $z2)\t// 33",
                "\t($x2 $y2 $z2)\t// 34",
                "\t($x3 $y2 $z2)\t// 35",
                " ",
                "\t($x0 $y0 $z3)\t// 36",
                "\t($x1 $y0 $z3)\t// 37",
                "\t($x2 $y0 $z3)\t// 38",
                "\t($x3 $y0 $z3)\t// 39",
                " ",
                "\t($x0 $y1 $z3)\t// 40",
                "\t($x1 $y1 $z3)\t// 41",
                "\t($x2 $y1 $z3)\t// 42",
                "\t($x3 $y1 $z3)\t// 43",
                " ",
                "\t($x0 $y2 $z3)\t// 44",
                "\t($x1 $y2 $z3)\t// 45",
                "\t($x2 $y2 $z3)\t// 46",
                "\t($x3 $y2 $z3)\t// 47",
            ");",
            " ",
            "//number of cells",
            "cx0\t40;",
            "cx1\t5;",
            " ",
            "cy0\t20;",
            "cy1\t40;",
            " ",
            "cz0\t30;",
            "cz1\t15;",
            " ",
            "//grading",
            "gx0\t-30;",
            "gx1\t30;",
            "gx2\t1;",
            " ",
            "gy0\t1;",
            "gy1\t2;",
            " ",
            "gz0\t-8;",
            "gz1\t8;",
            "gz2\t1;",
            " ",
            "blocks",
            "(",
                "\t//front",
                "\t//bottom front right //0",
                "\thex (0 1 5 4 12 13 17 16)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx0 $gy0 $gz0)",
                " ",
                "\t//bottom front centre //1",
                "\thex (12 13 17 16 24 25 29 28)",
                "\t($cx0 $cy0 $cz1)",
                "\tsimpleGrading ($gx0 $gy0 $gz2)",
                " ",
                "\t//bottom front left //2",
                "\thex (24 25 29 28 36 37 41 40)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx0 $gy0 $gz1)",
                " ",
                "\t//top front right //3",
                "\thex (4 5 9 8 16 17 21 20)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx0 $gy1 $gz0)",
                " ",
                "\t//top front centre //4",
                "\thex (16 17 21 20 28 29 33 32)",
                "\t($cx0 $cy1 $cz1)",
                "\tsimpleGrading ($gx0 $gy1 $gz2)",
                " ",
                "\t//top front left //5",
                "\thex (28 29 33 32 40 41 45 44)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx0 $gy1 $gz1)",
                " ",
                "\t//back",
                "\t//bottom back right //6",
                "\thex (2 3 7 6 14 15 19 18)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx1 $gy0 $gz0)",
                " ",
                "\t//bottom back centre //7",
                "\thex (14 15 19 18 26 27 31 30)",
                "\t($cx0 $cy0 $cz1)",
                "\tsimpleGrading ($gx1 $gy0 $gz2)",
                " ",
                "\t//bottom back left //8",
                "\thex (26 27 31 30 38 39 43 42)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx1 $gy0 $gz1)",
                " ",
                "\t//top back right //9",
                "\thex (6 7 11 10 18 19 23 22)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx1 $gy1 $gz0)",
                " ",
                "\t//top back centre //10",
                "\thex (18 19 23 22 30 31 35 34)",
                "\t($cx0 $cy1 $cz1)",
                "\tsimpleGrading ($gx1 $gy1 $gz2)",
                " ",
                "\t//top back left //11",
                "\thex (30 31 35 34 42 43 47 46)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx1 $gy1 $gz1)",
                " ",
                "\t//flap //12",
                "\t//flap bottom right",
                "\thex (1 2 6 5 13 14 18 17)",
                "\t($cx1 $cy0 $cz0)",
                "\tsimpleGrading ($gx2 $gy0 $gz0)",
                " ",
                "\t/*//flap bottom centre",
                "\thex (13 14 18 17 25 26 30 29)",
                "\t($cx1 $cy1 $cz1)",
                "\tsimpleGrading ($gx2 $gy0 $gz2)*/",
                " ",
                "\t//flap bottom left //13",
                "\thex (25 26 30 29 37 38 42 41)",
                "\t($cx1 $cy0 $cz0)",
                "\tsimpleGrading ($gx2 $gy0 $gz1)",
                " ",
                "\t//flap top right //14",
                "\thex (5 6 10 9 17 18 22 21)",
                "\t($cx1 $cy1 $cz0)",
                "\tsimpleGrading ($gx2 $gy1 $gz0)",
                " ",
                "\t//flap top centre //15",
                "\thex (17 18 22 21 29 30 34 33)",
                "\t($cx1 $cy1 $cz1)",
                "\tsimpleGrading ($gx2 $gy1 $gz2)",
                " ",
                "\t//flap top left //16",
                "\thex (29 30 34 33 41 42 46 45)",
                "\t($cx1 $cy1 $cz0)",
                "\tsimpleGrading ($gx2 $gy1 $gz1)",
            ");",
            " ",
            "boundary",
            "(",
                "\tinlet",
                "\t{",
                    "\t\ttype\tpatch;",
                    "\t\tfaces",
                    "\t\t(",
                        "\t\t\t(0 4 16 12)",
                        "\t\t\t(12 16 28 24)",
                        "\t\t\t(24 28 40 36)",
                        " ",
                        "\t\t\t(4 8 20 16)",
                        "\t\t\t(16 20 32 28)",
                        "\t\t\t(28 32 44 40)",
                    "\t\t);",
                "\t}",
                "\toutlet",
                "\t{",
                    "\t\ttype\tpatch;",
                    "\t\tfaces",
                    "\t\t(",
                        "\t\t\t(3 7 19 15)",
                        "\t\t\t(15 19 31 27)",
                        "\t\t\t(27 31 43 39)",
                        " ",
                        "\t\t\t(7 11 23 19)",
                        "\t\t\t(19 23 35 31)",
                        "\t\t\t(31 35 47 43)",
                    "\t\t);",
                "\t}",
                "\tflap",
                "\t{",
                    "\t\ttype\twall;",
                    "\t\tfaces",
                    "\t\t(",
                        "\t\t\t(13 17 29 25)",
                        "\t\t\t(14 18 17 13)",
                        "\t\t\t(26 30 18 14)",
                        "\t\t\t(25 29 30 26)",
                        " ",
                        "\t\t\t(17 18 30 29)",
                    "\t\t);",
                "\t}",
                "\tatmosphere",
                "\t{",
                    "\t\ttype\twall;",
                    "\t\tfaces",
                    "\t\t(",
                        "\t\t\t(8 9 21 20)",
                        "\t\t\t(20 21 33 32)",
                        "\t\t\t(32 33 45 44)",
                        " ",
                        "\t\t\t(9 10 22 21)",
                        "\t\t\t(21 22 34 33)",
                        "\t\t\t(33 34 46 45)",
                        " ",
                        "\t\t\t(10 11 23 22)",
                        "\t\t\t(22 23 35 34)",
                        "\t\t\t(34 35 47 46)",
                    "\t\t);",
                "\t}",
                "\tlowerWall",
                "\t{",
                    "\t\ttype\twall;",
                    "\t\tfaces",
                    "\t\t(",
                        "\t\t\t(0 1 13 12)",
                        "\t\t\t(12 13 25 24)",
                        "\t\t\t(24 25 37 36)",
                        " ",
                        "\t\t\t(1 2 14 13)",
                        "\t\t\t(25 26 38 37)",
                        " ",
                        "\t\t\t(2 3 15 14)",
                        "\t\t\t(14 15 27 26)",
                        "\t\t\t(26 27 39 38)",
                    "\t\t);",
                "\t}",
                "\tfrontAndBack",
                "\t{",
                    "\t\ttype\twall;",
                    "\t\tfaces",
                    "\t\t(",
                        "\t\t\t(0 1 5 4)",
                        "\t\t\t(1 2 6 5)",
                        "\t\t\t(2 3 7 6)",
                        " ",
                        "\t\t\t(4 5 9 8)",
                        "\t\t\t(5 6 10 9)",
                        "\t\t\t(6 7 11 10)",
                        " ",
                        "\t\t\t(36 37 41 40)",
                        "\t\t\t(37 38 42 41)",
                        "\t\t\t(38 39 43 42)",
                        " ",
                        "\t\t\t(40 41 45 44)",
                        "\t\t\t(41 42 46 45)",
                        "\t\t\t(42 43 47 46)",
                    "\t\t);",
                "\t}",
            ");",
        };
        
        printContent(outfile,content);
    }
}

void writeSolidPreciceDict(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        string filename = "preciceDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/solid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);
    
        vector<string> content = {
                "preciceConfig\t\"../precice-config.xml\";",
                " ",
                "participant Solid;",
                " ",
                "modules (FSI);",
                " ",
                "interfaces",
                "{",
                    "\tInterface1",
                    "\t{",
                        "\t\tmesh\tSolid-Mesh;",
                        "\t\tpatches\t(flap);",
                        "\t\tlocations\tfaceCenters;",
                        " ",
                        "\t\treadData",
                        "\t\t(",
                        "\t\t\tForce",
                        "\t\t);",
                        " ",
                        "\t\twriteData",
                        "\t\t(",
                        "\t\t\tDisplacement",
                        "\t\t);",
                    "\t};",
                "};",
                " ",
                "FSI",
                "{",
                    "\tsolverType\tsolid;",
                    " ",
                    "\tnameCellDisplacement\tD; //name of cell displacement field",
                    " ",
                    "\tnamePointDisplacement\tunused; //the solidDisplacementFoam does not have a point displacement field so we use the special name <unused>, which tells the adapter that it is not used",
                    " ",
                    "\tnameForce\tsolidForce; //name of the force field on the solid",
                "};"
            };
        
        printContent(outfile,content);
    }
}

void writeSolidControlDict(const string workDir, const vector<double>U, const vector<string>timestep, const vector<string>dirNames, string writeCompression){
    //using U as a counter
    for(int i = 0; i<U.size();i++){
        
        ostringstream strsCD;
        strsCD.clear();
        strsCD << "deltaT\t" << timestep.at(i) << ";";
        string valueCD;
        valueCD = strsCD.str();

        string filename = "controlDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/solid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "libs\t(\"libsolidDisplacementFoamForce.so\");",
            " ",
            "application\tsolidDisplacementFoam;",
            " ",
            "startFrom\tstartTime;",
            " ",
            "startTime\t0;",
            " ",
            "stopAt\tendTime;",
            " ",
            "endTime\t10;",
            " ",
            valueCD,
            " ",
            "writeControl\tadjustableRunTime;",
            " ",
            "writeInterval\t1e-02;",
            " ",
            "purgeWrite\t0;",
            " ",
            "writeFormat\tascii;",
            " ",
            "writePrecision\t12;",
            " ",
            "writeCompression\t" + writeCompression + ";",
            " ",
            "timeFormat\tgeneral;",
            " ",
            "timePrecision\t6;",
            " ",
            "functions",
            "{",
                "\tpreCICE_Adapter",
                "\t{",
                    "\t\ttype\tpreciceAdapterFunctionObject;",
                    "\t\tlibs\t(\"libpreciceAdapterFunctionObject.so\");",
                "\t}",
            "}"
        };
        
        printContent(outfile,content);
    }
}

void writeSolidFvSolution(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        string filename = "fvSolution";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/solid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "solvers",
            "{",
                "\t\"(D|T)\"",
                "\t{",
                    "\t\tsolver\tGAMG;",
                    "\t\ttolerance\t1e-09;",
                    "\t\trelTol\t0.001;",
                    "\t\tminIter\t3;",
                    "\t\tsmoother\tDIC;",
                    "\t\tnPreSweeps\t1;",
                    "\t\tnPostSweeps\t2;",
                    "\t\tnFinestSweeps\t2;",
                    "\t\tscaleCorrection\ttrue;",
                    "\t\tdirectSolveCoarsestLevel\tfalse;",
                    "\t\tcacheAgglomeration\ton; //don't use with dynamic mesh refinement",
                    "\t\tnCellsInCoarsestLevel 100; //sqrt(no. of cells)",
                    "\t\tagglomerator\tfaceAreaPair;",
                    "\t\tmergeLevels\t1;",
                "\t}",
            "}",
            " ",
            "stressAnalysis",
            "{",
                "\tcompactNormalStress\tyes;",
                "\tnCorrectors\t300; //Note: The accuracy of the solution can be significantly improved by increasing this number of iterations, but this will impact the runtime.",
                "\tD\t1e-06;",
            "}",
        };
        
        printContent(outfile,content);
    }
}

void writeSolidFvSchemes(const string workDir, const vector<double>U, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        string filename = "fvSchemes";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/solid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "d2dt2Schemes",
            "{",
                "\tdefault\tEuler;",
            "}",
            "ddtSchemes",
            "{",
                "\tdefault\tbackward;",
            "}",
            " ",
            "gradSchemes",
            "{",
                "\tdefault\tfourth;",
            "}",
            " ",
            "divSchemes",
            "{",
                "\tdefault\tnone;",
                "\tdiv(sigmaD)\tGauss cubic;",
            "}",
            " ",
            "laplacianSchemes",
            "{",
                "\tdefault\tnone;",
                "\tlaplacian(DD,D)\tGauss linear corrected;",
                "\tlaplacian(DT,T)\tGauss linear corrected;",
            "}",
            " ",
            "interpolationSchemes",
            "{",
                "\tdefault\tlinear;",
            "}",
            " ",
            "snGradSchemes",
            "{",
                "\tdefault\tnone;",
            "}"
        };
        
        printContent(outfile,content);
    }
}

void writeSolidDecomposeParDict(const string workDir, const vector<double>U, const vector<int>solidCores, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsSDPD;
        strsSDPD.clear();
        strsSDPD << solidCores.at(i);
        string valueSDPD;
        valueSDPD = strsSDPD.str();

        string filename = "decomposeParDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/solid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "numberOfSubdomains\t" + valueSDPD + ";",
            " ",
            "method\tscotch;"
        };
        
        printContent(outfile,content);
    }
}

void writeSolidBlockMeshDict(const string workDir, const vector<double>U, const vector<double>geometryProperties, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsBMD;
        strsBMD.clear();

        strsBMD << "x0\t" << -1*1e-3*geometryProperties.at(0)*geometryProperties.at(7)*0.5 << ";" << endl << "x1\t" << -1*1e-3*geometryProperties.at(3)*0.5 << ";" << endl << "x2\t" << 1e-3*geometryProperties.at(3)*0.5 << ";" << endl << "x3\t" << 1e-3*geometryProperties.at(0)*geometryProperties.at(7)*0.5 << ";" << endl << endl << "y0\t0" << ";" << endl << "y1\t" << 1e-3*geometryProperties.at(4) << ";" << endl << "y2\t" << 1e-3*geometryProperties.at(6) << ";" << endl << endl << "z0\t" << -1*1e-3*geometryProperties.at(2)*0.5 << ";" << endl << "z1\t" << -1*1e-3*geometryProperties.at(5)*0.5 << ";" << endl << "z2\t" << 1e-3*geometryProperties.at(5)*0.5 << ";" << endl << "z3\t" << 1e-3*geometryProperties.at(2)*0.5 << ";";

        string valueBMD;
        valueBMD = strsBMD.str();

        string filename = "blockMeshDict";
        string location = "system";    
        string filepath = "./" + workDir + "/." + location + "/" + dirNames.at(i) + "/solid/" + filename;
        string className = "dictionary";
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        printHeader(outfile,filename,location,className);

        vector<string> content = {
            "//flap base centred at 0,0,0",
            " ",
            valueBMD,
            " ",
            "vertices",
            "(",
                "\t($x0 $y0 $z0)\t// 0",
                "\t($x1 $y0 $z0)\t// 1",
                "\t($x2 $y0 $z0)\t// 2",
                "\t($x3 $y0 $z0)\t// 3",
                " ",
                "\t($x0 $y1 $z0)\t// 4",
                "\t($x1 $y1 $z0)\t// 5",
                "\t($x2 $y1 $z0)\t// 6",
                "\t($x3 $y1 $z0)\t// 7",
                " ",
                "\t($x0 $y2 $z0)\t// 8",
                "\t($x1 $y2 $z0)\t// 9",
                "\t($x2 $y2 $z0)\t// 10",
                "\t($x3 $y2 $z0)\t// 11",
                " ",
                "\t($x0 $y0 $z1)\t// 12",
                "\t($x1 $y0 $z1)\t// 13",
                "\t($x2 $y0 $z1)\t// 14",
                "\t($x3 $y0 $z1)\t// 15",
                " ",
                "\t($x0 $y1 $z1)\t// 16",
                "\t($x1 $y1 $z1)\t// 17",
                "\t($x2 $y1 $z1)\t// 18",
                "\t($x3 $y1 $z1)\t// 19",
                " ",
                "\t($x0 $y2 $z1)\t// 20",
                "\t($x1 $y2 $z1)\t// 21",
                "\t($x2 $y2 $z1)\t// 22",
                "\t($x3 $y2 $z1)\t// 23",
                " ",
                "\t($x0 $y0 $z2)\t// 24",
                "\t($x1 $y0 $z2)\t// 25",
                "\t($x2 $y0 $z2)\t// 26",
                "\t($x3 $y0 $z2)\t// 27",
                " ",
                "\t($x0 $y1 $z2)\t// 28",
                "\t($x1 $y1 $z2)\t// 29",
                "\t($x2 $y1 $z2)\t// 30",
                "\t($x3 $y1 $z2)\t// 31",
                " ",
                "\t($x0 $y2 $z2)\t// 32",
                "\t($x1 $y2 $z2)\t// 33",
                "\t($x2 $y2 $z2)\t// 34",
                "\t($x3 $y2 $z2)\t// 35",
                " ",
                "\t($x0 $y0 $z3)\t// 36",
                "\t($x1 $y0 $z3)\t// 37",
                "\t($x2 $y0 $z3)\t// 38",
                "\t($x3 $y0 $z3)\t// 39",
                " ",
                "\t($x0 $y1 $z3)\t// 40",
                "\t($x1 $y1 $z3)\t// 41",
                "\t($x2 $y1 $z3)\t// 42",
                "\t($x3 $y1 $z3)\t// 43",
                " ",
                "\t($x0 $y2 $z3)\t// 44",
                "\t($x1 $y2 $z3)\t// 45",
                "\t($x2 $y2 $z3)\t// 46",
                "\t($x3 $y2 $z3)\t// 47",
            ");",
            " ",
            "//number of cells",
            "cx0\t40;",
            "cx1\t5;",
            " ",
            "cy0\t20;",
            "cy1\t40;",
            " ",
            "cz0\t30;",
            "cz1\t15;",
            " ",
            "//grading",
            "gx0\t-30;",
            "gx1\t30;",
            "gx2\t1;",
            " ",
            "gy0\t1;",
            "gy1\t2;",
            " ",
            "gz0\t-8;",
            "gz1\t8;",
            "gz2\t1;",
            " ",
            "blocks",
            "(",
                "\t/*",
                "\t//front",
                "\t//bottom front right //0",
                "\thex (0 1 5 4 12 13 17 16)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx0 $gy0 $gz0)",
                " ",
                "\t//bottom front centre //1",
                "\thex (12 13 17 16 24 25 29 28)",
                "\t($cx0 $cy0 $cz1)",
                "\tsimpleGrading ($gx0 $gy0 $gz2)",
                " ",
                "\t//bottom front left //2",
                "\thex (24 25 29 28 36 37 41 40)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx0 $gy0 $gz1)",
                " ",
                "\t//top front right //3",
                "\thex (4 5 9 8 16 17 21 20)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx0 $gy1 $gz0)",
                " ",
                "\t//top front centre //4",
                "\thex (16 17 21 20 28 29 33 32)",
                "\t($cx0 $cy1 $cz1)",
                "\tsimpleGrading ($gx0 $gy1 $gz2)",
                " ",
                "\t//top front left //5",
                "\thex (28 29 33 32 40 41 45 44)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx0 $gy1 $gz1)",
                " ",
                "\t//back",
                "\t//bottom back right //6",
                "\thex (2 3 7 6 14 15 19 18)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx1 $gy0 $gz0)",
                " ",
                "\t//bottom back centre //7",
                "\thex (14 15 19 18 26 27 31 30)",
                "\t($cx0 $cy0 $cz1)",
                "\tsimpleGrading ($gx1 $gy0 $gz2)",
                " ",
                "\t//bottom back left //8",
                "\thex (26 27 31 30 38 39 43 42)",
                "\t($cx0 $cy0 $cz0)",
                "\tsimpleGrading ($gx1 $gy0 $gz1)",
                " ",
                "\t//top back right //9",
                "\thex (6 7 11 10 18 19 23 22)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx1 $gy1 $gz0)",
                " ",
                "\t//top back centre //10",
                "\thex (18 19 23 22 30 31 35 34)",
                "\t($cx0 $cy1 $cz1)",
                "\tsimpleGrading ($gx1 $gy1 $gz2)",
                " ",
                "\t//top back left //11",
                "\thex (30 31 35 34 42 43 47 46)",
                "\t($cx0 $cy1 $cz0)",
                "\tsimpleGrading ($gx1 $gy1 $gz1)",
                " ",
                "\t//flap //12",
                "\t//flap bottom right",
                "\thex (1 2 6 5 13 14 18 17)",
                "\t($cx1 $cy0 $cz0)",
                "\tsimpleGrading ($gx2 $gy0 $gz0)",
                "\t*/",
                " ",
                "\t//flap bottom centre",
                "\thex (13 14 18 17 25 26 30 29)",
                "\t($cx1 $cy1 $cz1)",
                "\tsimpleGrading ($gx2 $gy0 $gz2)",
                " ",
                "\t/*"
                "\t//flap bottom left //13",
                "\thex (25 26 30 29 37 38 42 41)",
                "\t($cx1 $cy0 $cz0)",
                "\tsimpleGrading ($gx2 $gy0 $gz1)",
                " ",
                "\t//flap top right //14",
                "\thex (5 6 10 9 17 18 22 21)",
                "\t($cx1 $cy1 $cz0)",
                "\tsimpleGrading ($gx2 $gy1 $gz0)",
                " ",
                "\t//flap top centre //15",
                "\thex (17 18 22 21 29 30 34 33)",
                "\t($cx1 $cy1 $cz1)",
                "\tsimpleGrading ($gx2 $gy1 $gz2)",
                " ",
                "\t//flap top left //16",
                "\thex (29 30 34 33 41 42 46 45)",
                "\t($cx1 $cy1 $cz0)",
                "\tsimpleGrading ($gx2 $gy1 $gz1)",
                "\t*/",
            ");",
            " ",
            "patches",
            "(",
                "\tpatch\tflap",
                "\t(",
                    "\t\t(13 17 29 25)",
                    "\t\t(14 18 17 13)",
                    "\t\t(26 30 18 14)",
                    "\t\t(25 29 30 26)",
                    " ",
                    "\t\t(17 18 30 29)",
                "\t)",
                "\tpatch\tbottom",
                "\t(",
                    "\t\t(13 14 26 25)",
                "\t)",
            ");",           
        };
        
        printContent(outfile,content);
    }
}

void writePreciceConfig(const string workDir, const vector<double>U, const vector<string>timestep, const vector<string>dirNames, double maxTime, const vector<double>geometryProperties){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsPC1;
        ostringstream strsPC2;
        strsPC1.clear();
        strsPC2.clear();
        strsPC1 << maxTime;

        //create 11 watchpoints along flap center
        for(int i=0;i<11;i++){
            strsPC2 << "\t\t<watch-point mesh=\"Solid-Mesh\" name=\"Flap-Line"<< i <<"\" coordinate=\"0;" << geometryProperties.at(4)/1000*0.1*i << ";0\"/>" << endl;
        }
        string valuePC1;
        string valuePC2;
        valuePC1 = strsPC1.str();
        valuePC2 = strsPC2.str();

        string filename = "precice-config.xml";
        string filepath = "./" + workDir + "/.precice-config/" + dirNames.at(i) + "/" + filename;
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        vector<string> content = {
            "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>",
            "<precice-configuration>",
                "\t<log enabled=\"true\">",
                    "\t\t<sink",
                        "\t\t\tfilter=\"%Severity% > debug and %Rank% = 0\"",
                        "\t\t\tformat=\"---[precice] %ColorizedSeverity% %Message%\"",
                        "\t\t\tenabled=\"true\"/>",
                "\t</log>",
                " ",
                "\t<data:vector name=\"Force\"/>",
                "\t<data:vector name=\"Displacement\"/>",
                " ",
                "\t<mesh name=\"Fluid-Mesh\" dimensions=\"3\">",
                    "\t\t<use-data name=\"Force\"/>",
                    "\t\t<use-data name=\"Displacement\"/>",
                "\t</mesh>",
                " ",
                "\t<mesh name=\"Solid-Mesh\" dimensions=\"3\">",
                    "\t\t<use-data name=\"Force\"/>",
                    "\t\t<use-data name=\"Displacement\"/>",
                "\t</mesh>",
                " ",
                "\t<participant name=\"Fluid\">",
                    "\t\t<read-data name=\"Displacement\" mesh=\"Fluid-Mesh\"/>",
                    "\t\t<write-data name=\"Force\" mesh=\"Fluid-Mesh\"/>",
                    "\t\t<provide-mesh name=\"Fluid-Mesh\"/>",
                    "\t\t<receive-mesh name=\"Solid-Mesh\" from=\"Solid\"/>",
                    "\t\t<mapping:nearest-neighbor direction=\"read\" from=\"Solid-Mesh\" to=\"Fluid-Mesh\" constraint=\"consistent\"/>",
                    "\t\t<mapping:nearest-neighbor direction=\"write\" from=\"Fluid-Mesh\" to=\"Solid-Mesh\" constraint=\"conservative\"/>",
                "\t</participant>",
                " ",
                "\t<participant name=\"Solid\">",
                    "\t\t<read-data name=\"Force\" mesh=\"Solid-Mesh\"/>",
                    "\t\t<write-data name=\"Displacement\" mesh=\"Solid-Mesh\"/>",
                    "\t\t<provide-mesh name=\"Solid-Mesh\"/>",
                    valuePC2,
                "\t</participant>",
                " ",
                "\t<m2n:sockets port=\"0\" acceptor=\"Solid\" connector=\"Fluid\" exchange-directory=\"../\" network=\"lo\" enforce-gather-scatter=\"0\" use-two-level-initialization=\"0\"/>",
                " ",
                "\t<coupling-scheme:parallel-implicit>",
                    "\t\t<participants first=\"Fluid\" second=\"Solid\"/>",
                    "\t\t<max-time value=\"" + valuePC1 + "\"/>",
                    "\t\t<time-window-size value=\"" + timestep.at(i) + "\"/>",
                    "\t\t<exchange data=\"Force\" mesh=\"Solid-Mesh\" from=\"Fluid\" to=\"Solid\"/>",
                    "\t\t<exchange data=\"Displacement\" mesh=\"Solid-Mesh\" from=\"Solid\" to=\"Fluid\"/>",
                    "\t\t<max-iterations value=\"100\"/>",
                    "\t\t<absolute-or-relative-convergence-measure abs-limit=\"1e-04\" rel-limit=\"5e-3\" data=\"Displacement\" mesh=\"Solid-Mesh\" strict=\"1\"/>",
                    "\t\t<absolute-or-relative-convergence-measure abs-limit=\"1e-04\" rel-limit=\"5e-3\" data=\"Force\" mesh=\"Solid-Mesh\" strict=\"1\"/>",
                    " ",
                    "\t\t<acceleration:IQN-ILS>",
                        "\t\t\t<data name=\"Displacement\" mesh=\"Solid-Mesh\"/>",
                        "\t\t\t<data name=\"Force\" mesh=\"Solid-Mesh\"/>",
                        "\t\t\t<preconditioner type=\"residual-sum\"/>",
                        "\t\t\t<filter type=\"QR2\" limit=\"1e-3\"/>",
                        "\t\t\t<initial-relaxation value=\"0.1\"/>",
                        "\t\t\t<max-used-iterations value=\"15\"/>",
                        "\t\t\t<time-windows-reused value=\"5\"/>",
                    "\t\t</acceleration:IQN-ILS>",
                "\t</coupling-scheme:parallel-implicit>",
            "</precice-configuration>"
        };
        
        for(const string &line : content){
            outfile << line << endl;
        }
    }
}

void writeInstallPreciceAndOF(){

    string filename = "installPreciceAndOF.sh";
    string filepath = filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    vector<string> content = {
        "#!/bin/bash -i",
        "#set installation directory",
        "installDir=~",
        " ",
        "#go to installation directory",
        " ",
        "cd $installDir",
        " ",
        "sudo apt-get -y update && sudo apt-get -y upgrade",
        "#install latest precice",
        "git clone https://github.com/precice/precice.git",
        "cd precice",
        "git checkout develop",
        "sudo apt -y install build-essential cmake libeigen3-dev libxml2-dev libboost-all-dev petsc-dev mercurial gnuplot pipx python3-gi python3-numpy python3-gi-cairo libgirepository1.0-dev gcc libcairo2-dev pkg-config python3-dev",
        "mkdir build",
        "cd build",
        "cmake -DCMAKE_BUILD_TYPE=Release ..",
        "make -j $(nproc)",
        "ctest ",
        "sudo make install",
        " ",
        "#these should not return empty",
        "pkg-config --cflags libprecice",
        "pkg-config --libs   libprecice",
        " ",
        "#install openfoam",
        "cd $installDir",
        "wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | sudo bash",
        "sudo apt -y install openfoam2312-dev",
        "echo \"source /usr/lib/openfoam/openfoam2312/etc/bashrc\" >> ~/.bashrc",
        " ",
        "#install latest openfoam adapter",
        "cd $installDir",
        "git clone https://github.com/precice/openfoam-adapter",
        "cd openfoam-adapter",
        "git checkout develop",
        "source ~/.bashrc",
        "./Allwmake -j",
        " ",
        "#download precice tutorials to build boundary condition",
        "cd $installDir/precice",
        "git clone https://github.com/precice/tutorials.git",
        "cd tutorials/perpendicular-flap/solid-openfoam/solidDisplacementFoamForce",
        "wmake libso",
        " ",
        "#install precice config viewer",
        "#use as: precice-config-visualizer precice-config.xml | dot -Tpdf > graph.pdf",
        "pipx install precice-config-visualizer",
        "pipx ensurepath"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writeBuildDirs(const string workDir, const vector<string>dirNames){

    string filename = "buildDirs.sh";
    string filepath = filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    string listOfDirs;
    string firstDir = dirNames.at(0);

    for(int i=1; i<dirNames.size();i++){
        if(i<dirNames.size()-1){
            listOfDirs += dirNames.at(i) + ",";
        }
        else{
            listOfDirs += dirNames.at(i);
        }
    }

    string allDirs = firstDir + "," + listOfDirs;

    vector<string> content = {
        "#!/bin/bash -i",
        "#build top directory",
        "mkdir " + workDir,
        "cd " + workDir,
        " ",
        "mkdir {cases,.0,.constant,.system,.precice-config,.scripts}",
        " ",
        "mkdir cases/" + firstDir,
        "mkdir cases/" + firstDir + "/fluid",
        "mkdir cases/" + firstDir + "/fluid/{0,constant,system}",
        "cp -r cases/" + firstDir + "/fluid cases/" + firstDir + "/solid",
        "touch cases/" + firstDir + "/{fluid/fluid.foam,solid/solid.foam}",
        "for i in {" + listOfDirs + "}; do cp -r cases/" + firstDir + " cases/\"$i\"; done",
        " ",
        "mkdir .0/{" + allDirs + "}",
        " ",
        "mkdir .constant/{" + allDirs + "}",
        " ",
        "mkdir .system/"+ firstDir,
        "mkdir .system/" + firstDir + "/fluid",
        "mkdir .system/" + firstDir + "/solid",
        "for i in {" + listOfDirs + "}; do cp -r .system/" + firstDir + " .system/\"$i\"; done",
        " ",
        "mkdir .precice-config/{" + allDirs + "}",
        " ",
        "mkdir .scripts/{" + allDirs + "}"

    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writeAllrun(const string executableName){

    string filename = "allrun.sh";
    string filepath = filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    vector<string> content = {
        "#!/bin/bash -i",
        "./installPreciceAndOF.sh",
        "./buildDirs.sh",
        executableName + " build",
        "./copyFiles.sh"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writeCopyFiles(const string workDir, const vector<string>dirNames){

    string filename = "copyFiles.sh";
    string filepath = filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    string allDirs;

    for(int i=0; i<dirNames.size();i++){
        if(i<dirNames.size()-1){
            allDirs += dirNames.at(i) + ",";
        }
        else{
            allDirs += dirNames.at(i);
        }
    }

    vector<string> content = {
        "#!/bin/bash -i",
        "cd " + workDir,
        "for i in {" + allDirs + "}; do cp -r .0/\"$i\"/{k,nut,omega,p,phi,pointDisplacement,U} cases/\"$i\"/fluid/0; done",
        "for i in {" + allDirs + "}; do cp -r .0/\"$i\"/{D,solidForce} cases/\"$i\"/solid/0; done",
        "for i in {" + allDirs + "}; do cp -r .constant/\"$i\"/{dynamicMeshDict,transportProperties,turbulenceProperties} cases/\"$i\"/fluid/constant; done",
        "for i in {" + allDirs + "}; do cp -r .constant/\"$i\"/{mechanicalProperties,thermalProperties} cases/\"$i\"/solid/constant; done",
        "for i in {" + allDirs + "}; do cp -r .system/\"$i\"/fluid/{blockMeshDict,controlDict,decomposeParDict,fvSchemes,fvSolution,preciceDict} cases/\"$i\"/fluid/system; done",
        "for i in {" + allDirs + "}; do cp -r .system/\"$i\"/solid/{blockMeshDict,controlDict,decomposeParDict,fvSchemes,fvSolution,preciceDict} cases/\"$i\"/solid/system; done",
        "for i in {" + allDirs + "}; do cp -r .precice-config/\"$i\"/precice-config.xml cases/\"$i\"; done",
        "for i in {" + allDirs + "}; do cp -r .scripts/{clean.sh,runFluidCluster.sh,plotDisplacement.sh} cases/\"$i\"/fluid; done",
        "for i in {" + allDirs + "}; do cp -r .scripts/{clean.sh,plotDisplacement.sh} cases/\"$i\"/solid; done",
        "for i in {" + allDirs + "}; do cp -r .scripts/\"$i\"/{runFluid.sh,serverJobFluid} cases/\"$i\"/fluid; done",
        "for i in {" + allDirs + "}; do cp -r .scripts/\"$i\"/{runSolid.sh,runSolidCluster.sh,serverJobSolid} cases/\"$i\"/solid; done"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writeClean(const string workDir){

    string filename = "clean.sh";
    string filepath = "./" + workDir + "/.scripts/" + filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    vector<string> content = {
        "#!/bin/bash -i",
        "shopt -s extglob",
        "rm -rv !(\"0\"|\"constant\"|\"system\"|\"fluid.foam\"|\"solid.foam\"|\"runFluid.sh\"|\"runSolid.sh\"|\"runFluidCluster.sh\"|\"runSolidCluster.sh\"|\"clean.sh\"|\"plotDisplacement.sh\"|\"serverJobFluid\"|\"serverJobSolid\")",
        "rm -r constant/polyMesh"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writeRunFluid(const string workDir, const vector<double>U, const vector<int>fluidCores, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsRF;
        strsRF.clear();
        strsRF << fluidCores.at(i);
        string valueRF;
        valueRF = strsRF.str();

        string filename = "runFluid.sh";
        string filepath = "./" + workDir + "/.scripts/"  + dirNames.at(i) + "/" + filename;
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        vector<string> content = {
            "#!/bin/bash -i",
            "blockMesh",
            "checkMesh | tee log.checkMesh",
            "decomposePar -force",
            "mpirun -np " + valueRF + " pimpleFoam -parallel | tee log.pimple"
        };
        
        for(const string &line : content){
            outfile << line << endl;
        }
    }
}

void writeRunSolid(const string workDir, const vector<double>U, const vector<int>solidCores, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsRS1;
        ostringstream strsRS2;
        strsRS1.clear();
        strsRS1 << solidCores.at(i) - 1;
        string valueRS1;
        string valueRS2;
        valueRS1 = strsRS1.str();
        strsRS2.clear();
        strsRS2 << solidCores.at(i);
        valueRS2 = strsRS2.str();

        string filename = "runSolid.sh";
        string filepath = "./" + workDir + "/.scripts/" + dirNames.at(i) + "/" + filename;
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        vector<string> content = {
            "#!/bin/bash -i",
            "blockMesh",
            "checkMesh | tee log.checkMesh",
            "decomposePar -force",
            "#correct wrong 0 dir files in proc* from badly decomposed D BC",
            "for i in {0.." + valueRS1 + "..1}",
            "do ",
                "\tcd ./processor\"$i\"/0",
                "\tsed -i \"/force           uniform (0 0 0);/,/gradient        uniform (0 0 0);/d\" D",
                "\tsed -i \"/forceField	solidForce;/,/value	uniform (0 0 0);/d\" D",
                "\tsed -i \"/type            solidDisplacementFoamForce;/a\\  \\t\\tforceField\\tsolidForce;\\n\\t\\tvalue\\tuniform (0 0 0);\" D",
                "\tcd ../../",
            "done",
            " ",
            "mpirun -np " + valueRS2 + " solidDisplacementFoam -parallel | tee log.solid"
        };
        
        for(const string &line : content){
            outfile << line << endl;
        }
    }
}

void writeRunFluidCluster(const string workDir){

    string filename = "runFluidCluster.sh";
    string filepath = "./" + workDir + "/.scripts/" + filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    vector<string> content = {
        "#!/bin/bash",
        "module load openfoam/gcc/9.1.0/v2212",
        "source /cluster/applications/openfoam/gcc/9.1.0/OpenFOAM-v2212/etc/bashrc",
        " ",
        "blockMesh",
        "checkMesh | tee log.checkMesh",
        "decomposePar -force",
        "qsub serverJobFluid"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writeRunSolidCluster(const string workDir, const vector<double>U, const vector<int>solidCores, const vector<string>dirNames){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsRSC;
        strsRSC.clear();
        strsRSC << solidCores.at(i) - 1;
        string valueRSC;
        valueRSC = strsRSC.str();

        string filename = "runSolidCluster.sh";
        string filepath = "./" + workDir + "/.scripts/" + dirNames.at(i) + "/" + filename;
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        vector<string> content = {
            "#!/bin/bash",
            "module load openfoam/gcc/9.1.0/v2212",
            "source /cluster/applications/openfoam/gcc/9.1.0/OpenFOAM-v2212/etc/bashrc",
            " ",
            "blockMesh",
            "checkMesh | tee log.checkMesh",
            "decomposePar -force",
            "#correct wrong 0 dir files in proc* from badly decomposed D BC",
            "for i in {0.." + valueRSC + "..1}",
            "do ",
                "\tcd ./processor\"$i\"/0",
                "\tsed -i \"/force           uniform (0 0 0);/,/gradient        uniform (0 0 0);/d\" D",
                "\tsed -i \"/forceField	solidForce;/,/value	uniform (0 0 0);/d\" D",
                "\tsed -i \"/type            solidDisplacementFoamForce;/a\\  \\t\\tforceField\\tsolidForce;\\n\\t\\tvalue\\tuniform (0 0 0);\" D",
                "\tcd ../../",
            "done",
            " ",
            "qsub serverJobSolid"
        };
        
        for(const string &line : content){
            outfile << line << endl;
        }
    }
}

void writeServerJobFluid(const string workDir, const vector<double>U, const vector<int>fluidCores, const vector<string>dirNames, const string mailPush){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsSJF;
        strsSJF.clear();
        strsSJF << fluidCores.at(i);
        string valueSJF;
        valueSJF = strsSJF.str();

        string filename = "serverJobFluid";
        string filepath = "./" + workDir + "/.scripts/" + dirNames.at(i) + "/" + filename;
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        vector<string> content = {
            "#!/bin/bash",
            "#PBS -l select=64:ncpus=1:mem=1000mb,place=free",
            "#PBS -l walltime=4:00:00",
            "#PBS -m abe",
            "#PBS -M " + mailPush,
            "#PBS -N imfdFlap-fluid",
            " ",
            "module load openfoam/gcc/9.1.0/v2212",
            "source /cluster/applications/openfoam/gcc/9.1.0/OpenFOAM-v2212/etc/bashrc"
            " ",
            "cd $PBS_O_WORKDIR",
            " ",
            "mpirun -np " + valueSJF + " pimpleFoam -parallel | tee log.pimple"
        };
        
        for(const string &line : content){
            outfile << line << endl;
        }
    }
}

void writeServerJobSolid(const string workDir, const vector<double>U, const vector<int>solidCores, const vector<string>dirNames, const string mailPush){
    //using U as a counter
    for(int i = 0; i<U.size();i++){

        ostringstream strsSJS;
        strsSJS.clear();
        strsSJS << solidCores.at(i);
        string valueSJS;
        valueSJS = strsSJS.str();

        string filename = "serverJobSolid";
        string filepath = "./" + workDir + "/.scripts/" + dirNames.at(i) + "/" + filename;
        ofstream outfile;
        outfile.open(filepath, ios::binary);

        vector<string> content = {
            //needs to rebuild solid BC on decompose
            "#!/bin/bash",
            "#PBS -l select=64:ncpus=1:mem=1000mb,place=free",
            "#PBS -l walltime=4:00:00",
            "#PBS -m abe",
            "#PBS -M " + mailPush,
            "#PBS -N imfdFlap-fluid",
            " ",
            "module load openfoam/gcc/9.1.0/v2212",
            "source /cluster/applications/openfoam/gcc/9.1.0/OpenFOAM-v2212/etc/bashrc"
            " ",
            "cd $PBS_O_WORKDIR",
            " ",
            "mpirun -np " + valueSJS + " solidDisplacementFoam -parallel | tee log.pimple"
        };
        
        for(const string &line : content){
            outfile << line << endl;
        }
    }
}

void writePlotDisplacement(const string workDir){

    string filename = "plotDisplacement.sh";
    string filepath = "./" + workDir + "/.scripts/" + filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    vector<string> content = {
        "#!/bin/sh",
        "gnuplot -p << EOF",
        "set grid",
        "set title 'x-displacement of the flap tip'",
        "set xlabel 'time [s]'",
        "set ylabel 'x-displacement [m]'",
        "#set term pngcairo enhanced size 900,654",
        "#set output \"all-watchpoints.png\"",
        "set term qt",
        "plot \"precice-Solid-watchpoint-Flap-Tip.log\" using 1:8 with lines title \"watchpoint-X\", \"precice-Solid-watchpoint-Flap-Tip.log\" using 1:9 with lines title \"watchpoint-Y\", \"precice-Solid-watchpoint-Flap-Tip.log\" using 1:10 with lines title \"watchpoint-Z\"",
        "while (1) {",
		"replot",
		"pause 5}",
        "EOF"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writeInputParameters(const string workDir){

    string filename = "inputParameters";
    string filepath = filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    vector<string> content = {
        "//This is the default inputParameters file generated because no file with the same name was found in the same directory as the executable.",
        "//lists volFlowRate, dirNames, timestep, fluidCores, solidCores can be changed in length to build more or fewer directories but they must all be of equal size.",
        "//lists fluidProperties, solidProperties, geometricProperties may not be changed in length without rewriting the source code.",
        "//set string for HPC mailpush or leave blank after keyword:",
        "mailPush alrik-matti.luda@student.tu-freiberg.de",
        "//set compression flag for file output of openfoam runs:",
        "writeCompression off",
        "//set maximum simulation time:",
        "maxTime 10",
        "//work directory is set to:",
        "workDir " + workDir,
        "//volumentric flow rates in m3/h are set to:",
        "volFlowRate 10,20,30,40,50,60,70,80,90,100",
        "//directory names will be set to:",
        "dirNames 10-m3,20-m3,30-m3,40-m3,50-m3,60-m3,70-m3,80-m3,90-m3,100-m3",
        "//timesteps in s are set to:",
        "timestep 1e-2,5e-3,2.5e-3,2.5e-3,2.5e-3,2e-3,1e-3,1e-3,1e-3,1e-3",
        "//number of cores the fluid solver will run on:",
        "fluidCores 16,16,16,16,16,16,16,16,16,16",
        "//number of cores the solid solver will run on:",
        "solidCores 8,8,8,8,8,8,8,8,8,8",
        "//poperties of the simulated fluid using the format <density in kg/m3>,<dynamic viscosity in kg/(m*s)>",
        "fluidProperties 998.2,0.0010003",
        "//poperties of the simulated solid using the format <density in kg/m3>,<Young's modulus E in Pa>,<Possion's ratio>,<shear modulus G in Pa>",
        "solidProperties 1020,1.23e6,0.3,4.73e5",
        "//properties of the water channel and the flap using the format <channel length L mm>,<channel height H mm>,<channel width W mm>,<lap length streamwise l mm>,<flap height h mm>,<flap width w mm>,<height of water in channel Hw mm>,<fraction of channel length as length for blockMesh domain>",
        "geometricProperties 5000,405,309,5,100,20,250,0.2"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}

void writePlotSolverTime(const string workDir){

    string filename = "plotSolverTime.sh";
    string filepath = "./" + workDir + "/.scripts/" + filename;
    ofstream outfile;
    outfile.open(filepath, ios::binary);

    vector<string> content = {
        "#!/bin/sh",
        "gnuplot -p << EOF",

        "set grid",
        "set title 'Total and solver run-time'",
        "set xlabel 'Iteration'",
        "set ylabel 'time [s]'",
        "#set term pngcairo enhanced size 900,654",
        "#set output \"all-watchpoints.png\"",
        "set term qt",
        "plot \"< cat log.pimple | grep 'ExecutionTime = ' | cut -d' ' -f3 | tr -d ','\" title 'Execution Time Pimple' with lines,\\",
        "\"<cat log.pimple | grep 'ExecutionTime = ' | cut -d' ' -f8 | tr -d ','\" title 'Clock Time' with lines,\\",
        "\"<cat ../solid/log.solid | grep 'ExecutionTime = ' | cut -d' ' -f3 | tr -d ','\" title 'Execution Time Solid' with lines",

        "while (1) {",
        "replot",
        "pause -5",
        "}",
        "EOF"
    };
        
    for(const string &line : content){
        outfile << line << endl;
    }

}


vector<double> readToVariableDouble(string input){
    vector<double>doubleVector;
    while(input.size()>0){
        if(input.find(",") != string::npos){
            size_t pos = input.find_first_of(",");
            string subStr = input.substr(0,pos);
            input.erase(input.begin(),input.begin()+pos+1);
            doubleVector.push_back(stod(subStr));
            subStr.erase();
        }
        else{
            doubleVector.push_back(stod(input));
            input.erase();
        }
    }
    return doubleVector;
}

vector<string> readToVariableString(string input){
    vector<string>stringVector;
    while(input.size()>0){
            if(input.find(",") != string::npos){
                size_t pos = input.find_first_of(",");
                string subStr = input.substr(0,pos);
                input.erase(input.begin(),input.begin()+pos+1);
                stringVector.push_back(subStr);
                subStr.erase();
            }
            else{
                stringVector.push_back(input);
                input.erase();
            }
        }
    return stringVector;
}

vector<int> readToVariableInt(string input){
    vector<int>intVector;
    while(input.size()>0){
        if(input.find(",") != string::npos){
            size_t pos = input.find_first_of(",");
            string subStr = input.substr(0,pos);
            input.erase(input.begin(),input.begin()+pos+1);
            intVector.push_back(stod(subStr));
            subStr.erase();
        }
        else{
            intVector.push_back(stoi(input));
            input.erase();
        }
    }
    return intVector;
}



int main(int argc, char * argv[]){
    
    string executableName = argv[0];

    string lastEdit = "2024-04-14";

    string keyword;

    string workDir = "imfdFlap";

    if(argc>1){
        keyword = argv[1];
    }

    //only run this with keyword build or prepare passed as argument
    if(keyword=="build" || keyword == "prepare"){

        //read inputParameters from file
        //check if inputParameters exists as a file, else write a default file to be read back in
            if(fs::exists("inputParameters")){
                cout << "Parameters read from inputParameters file." << endl;
            }
            else{
                cout << "Generating default inputParameters file." << endl;
                writeInputParameters(workDir);

            }

        vector<string>requiredParameters={"mailPush","writeCompression","workDir","maxTime","volFlowRate","dirNames","timestep","fluidCores","solidCores","fluidProperties","solidProperties","geometricProperties"};
        string inputParameterBuffer[requiredParameters.size()];


        ifstream infile("inputParameters");
        string readLine;
        while(getline(infile,readLine)){
            if(readLine.find("//") != string::npos){
                cout << "Dropping commented line: " + readLine << endl;
            }
            else{
                for(int i=0; i<requiredParameters.size();i++){
                    if(readLine.find(requiredParameters.at(i)) != string::npos){
                        cout << "Found " << requiredParameters.at(i) << endl;
                        int delCount =  requiredParameters.at(i).size() + 1; //delete keyword and whitespace

                        readLine.erase(readLine.begin(), readLine.begin() + delCount);
                        inputParameterBuffer[i] = readLine;
                    }
                }

            }
        }

        //return values from inputParameterBuffer to variables
        string mailPush = inputParameterBuffer[0];

        string writeCompression = inputParameterBuffer[1];

        workDir = inputParameterBuffer[2];

        double maxTime = stod(inputParameterBuffer[3]);


        vector<double>volFlowRate = readToVariableDouble(inputParameterBuffer[4]);

        vector<string>dirNames = readToVariableString(inputParameterBuffer[5]);
        vector<string>timestep = readToVariableString(inputParameterBuffer[6]);

        vector<int>fluidCores = readToVariableInt(inputParameterBuffer[7]);
        vector<int>solidCores = readToVariableInt(inputParameterBuffer[8]);

        vector<double>fluidProperties = readToVariableDouble(inputParameterBuffer[9]);
        vector<double>solidProperties = readToVariableDouble(inputParameterBuffer[10]);
        vector<double>geometryProperties = readToVariableDouble(inputParameterBuffer[11]);

        //check vectors are the same size or defined size
        if(volFlowRate.size() == dirNames.size() && volFlowRate.size() == timestep.size() && volFlowRate.size() == fluidCores.size() && volFlowRate.size() == solidCores.size() &&  fluidProperties.size() == 2 && solidProperties.size() == 4 && geometryProperties.size() == 8){
            cout << "Lists passed length check during "+ keyword +"." << endl << "Files will be written to " << volFlowRate.size() << " case directories." << endl;
        }
        else{
            cout << "Lists failed length check during " + keyword +". Either regenerate default inputParameters file or check edits to inputParameters. Make sure lists volFlowRate, dirNames, timestep, fluidCores, solidCores are all of equal size and lists fluidProperties, solidProperties, geometricProperties contain 2, 4, 8 elements respectively." << endl;
            cout << "volFlowRate contains " << volFlowRate.size() << " elements." << endl;
            cout << "dirNames contains " << dirNames.size() << " elements." << endl;
            cout << "timestep contains " << timestep.size() << " elements." << endl;
            cout << "fluidCores contains " << fluidCores.size() << " elements." << endl;
            cout << "solidCores contains " << solidCores.size() << " elements." << endl;
            cout << "fluidProperties contains " << fluidProperties.size() << " elements." << endl;
            cout << "solidProperties contains " << solidProperties.size() << " elements." << endl;
            cout << "geometryProperties contains " << geometryProperties.size() << " elements." << endl;
            return 0;

        }

        if(keyword=="build"){

            //vectors to make velocity and derived turbulence values available
            //no vectors for phi, p, pointDisplacement as these remain unaffected by a change in u

            vector<double>U;
            vector<double>omega;
            vector<double>k;
            vector<double>I;

            //calculate values and fill vectors U, omega and k
            for(int i=0;i<10;i++){
                //u in m/s from 10 m3/h -> 100 m3/h using A = W*Hw
                U.push_back((volFlowRate.at(i))/3600/(1e-6*geometryProperties.at(2)*geometryProperties.at(6)));
                //Reynolds number of channel with hydraulic diameter
                double Re = fluidProperties.at(0)*U.at(i)*geometryProperties.at(2)/(2*geometryProperties.at(6)+2*geometryProperties.at(2))/fluidProperties.at(1);
                //Turbulence intensity for fully developed pipe flow
                I.push_back(0.16*pow(Re,-0.125));
                //Turbulent kinetic energy
                k.push_back(1.5*I.at(i)*pow(U.at(i),2));
                //omega
                omega.push_back(pow(k.at(i),0.5)/(pow(0.09,0.25)*1e-3*geometryProperties.at(5)));

            }

            //O
            writeU(workDir,U,geometryProperties,dirNames);
            writeK(workDir,k,I,dirNames);
            writeOmega(workDir,omega,dirNames);
            writeP(workDir,U,dirNames);
            writeNut(workDir,U,dirNames);
            writePhi(workDir,U,dirNames);
            writePointDisplacement(workDir,U,dirNames);
            writeD(workDir,U,dirNames),
            writeSolidForce(workDir,U,dirNames);

            //constant
            writeDynamicMeshDict(workDir,U,dirNames);
            writeTransportProperties(workDir,U,fluidProperties,dirNames);
            writeTurbulenceProperties(workDir,U,dirNames);
            writeThermalProperties(workDir,U,dirNames);
            writeMechanicalProperties(workDir,U,solidProperties,dirNames);

            //system
            writeFluidControlDict(workDir,U,timestep,dirNames,writeCompression);
            writeFluidPreciceDict(workDir,U,fluidProperties,dirNames);
            writeFluidFvSolution(workDir,U,dirNames);
            writeFluidFvSchemes(workDir,U,dirNames);
            writeFluidDecomposeParDict(workDir,U,fluidCores,dirNames);
            writeFluidBlockMeshDict(workDir,U,geometryProperties,dirNames);

            writeSolidControlDict(workDir,U,timestep,dirNames,writeCompression);
            writeSolidPreciceDict(workDir,U,dirNames);
            writeSolidFvSolution(workDir,U,dirNames);
            writeSolidFvSchemes(workDir,U,dirNames);
            writeSolidDecomposeParDict(workDir,U,solidCores,dirNames);
            writeSolidBlockMeshDict(workDir,U,geometryProperties,dirNames);

            //precice-config
            writePreciceConfig(workDir,U,timestep,dirNames,maxTime,geometryProperties);

            //scripts
            writeClean(workDir);
            writeRunFluid(workDir,U,fluidCores,dirNames);
            writeRunSolid(workDir,U,solidCores,dirNames);
            writeRunFluidCluster(workDir);
            writeRunSolidCluster(workDir,U,solidCores,dirNames);
            writeServerJobFluid(workDir,U,fluidCores,dirNames,mailPush);
            writeServerJobSolid(workDir,U,solidCores,dirNames,mailPush);
            writePlotDisplacement(workDir);
        }
        else if(keyword=="prepare"){
            //check if new workDir name supplied
            if(keyword=="prepare" && argc==3){
                workDir = argv[2];
                cout << "workDir set to " << workDir << " via argument and will override inputParamenters if it already exists and specifies a different directory." << endl;
            }
            else{
                workDir = "imfdFlap";
            }

            //write on first run
            writeInstallPreciceAndOF();
            writeBuildDirs(workDir,dirNames);
            writeAllrun(executableName);
            writeCopyFiles(workDir,dirNames);

        }
        else{
            cout << "End of program." << endl;
            return 0;
        }
    }
    else if(keyword=="wcgw" || keyword == "help"){
        cout << "Welcome to the What Could Go Wrong? section of my file. If you turn up here, chances are you're miserably stuck so here's a list of issues I encountered that may help you unstick yourself.\n\nIf I haven't described the issue you're facing below, I may have wrongly passed it off as a neglible nuisance. Should this be the case, please contact me erstwhile at alrik-matti.luda@student.tu-freiberg.de until I can be bothered to finish setting up github.\n\n>>Inaccurate results:<<\nHaving adapted everything from the 2D perpendicular-flap tutorial found on the preCICE website (source: https://precice.org/tutorials-perpendicular-flap.html), I wasn't suprised to find almost every setting in the solid and fluid fvScheme set to some first order scheme. The only significant improvements in terms of reproducibility (compared to other solvers coupled to openfoam) were achieved by setting\n\ngradSchemes default -> fourth\n\ndivSchemes default -> Gauss linear\n\n>>Serial run works but parallel run fails with processorPolyPatch error:<<\nTake a look at the solution I described here:  https://www.cfd-online.com/Forums/openfoam-solving/255295-pimplefoam-aborts-during-parallel-run-processorpolypatch-error.html\n\n>>Serial run works but parallel run results in bad communication:<<\nCheck the precice-config and find out if you have a loopback device addressed by \"lo\" or if it has a different address.\n\n>>Managing a high Courant number:<<\nThe Courant number evolves over time as the mesh is moved and reflects the highest velocity in the smallest length scale. Slamming the fluid into the flap without a ramp may squish the cells beyond the point of no return. The solver will abort as soon as the cells crash into each other. Even if the Courant number declines after the initial bounce, you might now be facing a number significantly below 1 with no way to use adjustableTimeStep due to a lack of support for adjustableTimeStep by the openfoam-adapter. Tl;dr: use a ramp (default setting in the generated U BC).\n\n>>Coupling is successful but solidDisplacementFoam returns 0 on all calculations:<<\nThis may be caused by the custom force BC for solidDisplacementFoam not decomposing correctly. I implemented a workaround in the runSolid.sh script so I hope you will never encounter it. Essentially: copy the flap BC from 0/D to all proc*/0/D's.\n\n>>Can't resume simulations from latest time:<<\nprecice-config.xml needs to be set to a later time as well and the start time adjusted.\n\n>>How can I track the progress of my coupled simulation:<<\nBesides tailing the solver logs you can take a look at the convergence and iteration logs produced by precice or run the plot-displacement.sh script to visualize your progress so far.\n\n>>I want to run the turek-hron-fsi3 case but I don't have groovyBC installed:<<\nGlad you asked, here's an installation guide, which wasn't all that easy to find:\nDownload and install swak4foam, which is required by the turek-hron-fsi3 precice tutorial\ncd \"$HOME/OpenFOAM/$USER-$WM_PROJECT_VERSION\"\nhg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam swak4Foam\ncd swak4Foam\nhg update develop\necho \"Run ./AllwmakeAll -j -q -l inside the openfoam shell\"\nopenfoam2306\n\n>>I want to try deal.II:<<\nSure, here you go:\nsudo apt install libdeal.ii-dev libdeal.ii-doc cmake make g++\n#this is optional, build tutorials to test deal.ii\ncp -r /usr/share/doc/libdeal.ii-doc/examples/step-1 .\ncd step-1\ncmake .\nmake run\n#the precice-dealii-adapter\ngit clone https://github.com/precice/dealii-adapter.git && cd dealii-adapter\n#defaults to -DDIM=2 and builds a 2D solver called elasticity\ncmake . //cmake -DDIM=3 .\nmake\n#add the elasticity executables path to your environment or just copy it into the case directory\n\n>>I want to try FEniCs:<<\nSure, here you go:\nsudo apt install software-properties-common\nsudo add-apt-repository ppa:fenics-packages/fenics\nsudo apt update\nsudo apt install fenicsx\nsudo apt install python3-pip\npip3 install --user fenicsprecice\n\n>>I want to try FEniCs:<<\nSure, here you go:\n#change lines in perpendicular-flap tutorial case because ufl has been removed from fenics-packages/fenics\nchange in the solid.py:\nreplace: from ufl import nabla_div\nwith: from ufl_legacy import nabla_div\n\timport ufl_legacy as ufl\n\n>>I want to try Calculix:<<\nSure, but it will only work with preCICE v2 here you go:\nsudo apt install libarpack2-dev libspooles-dev libyaml-cpp-dev\nwget http://www.dhondt.de/ccx_2.20.src.tar.bz2\ntar xvjf ccx_2.20.src.tar.bz2\n#make sure you have the right package for your OS, this one is for jammy\nwget https://github.com/precice/calculix-adapter/releases/download/v2.20.0/calculix-precice2_2.20.0-1_amd64_jammy.deb\nsudo apt install ./calculix-precice2_2.20.0-1_amd64_jammy.deb\ncd CalculiX\nwget https://github.com/precice/calculix-adapter/archive/refs/heads/master.tar.gz\ntar -xzf master.tar.gz\n#fork to take care of a broken file\ncd ccx_2.20/src\nrm cubtri.f\n#retrieve replacement from https://calculix.discourse.group/t/compilation-error-with-gfortran10/1048/3\ncp /location/of/replacement/cubtri.f .\ncd ../..\n#end of fork\ncd calculix-adapter-master\nmake\nexport PATH=\"/home/ocelittle/CalculiX/calculix-adapter-master/bin:${PATH}\"\n\n>>The simulations occupy a lot of disk space:<<\nwriteCompression is off by default but I string to manage it is available in the main function.\n\n>>Where did you find your values and equations for your turbulence values:<<\nRight here and here:\nhttps://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras-k-omega-sst.html\n\nhttps://www.cfd-online.com/Wiki/Turbulence_intensity\nThese are the most major hurdles I had to take, but there are a lot more to trip over; unfortunately, I forgot to write them down. Add them yourself or message me.\n\n\n\nHappy coupling!\n\n-- Alrik Luda, Freiberg, 2024-04-10" << endl;
    }
    else{
        cout << "Self-installing OpenFOAM experiment: IMFD-Flap using preCICE, by Alrik Luda - rev. 0 - (" + lastEdit + ")\n\n>>About:<<\nThis file generates scripts to install the latest versions of OpenFOAM, OpenFOAM-Adapter, preCICE and anything else required to run them at the time of writing. It also builds a bunch of case directories (under cases) varying in their volumetric flow rate (from 10 m3 to 100 m3). To start a coupled simulation enter the solid and fluid directories of a case with two seperate terminals and execute the appropriate run script.\n\n>>Problems during installation:<<\nThis installation script was built on Ubuntu 22.04.4 LTS (Jammy Jellyfish) running on WSL on Windows 11. Verify the OpenFOAM repo (and others) are compatible or commence an edit in the installPreciceAndOF.sh file. If the installation fails, someone probably updated their version without telling anyone else, making it necessary to fall back to fixed versions compatible with one another.\n\nEverything is intended to run on linux, but this file can also be compiled on windows; just make sure the LF end-of-line flag is set or else all created files will be bogus.\n\n>>Just a heads-up:<<\nThis file contains a string mailPush, which is originally set to my e-mail address. Please keep this in mind before using the generated files to queue a HPC job - or we'll be in touch. Default walltime is set to 4h but will need to be adjusted depending on the available cores and hardware.\n\n>>Arguments:<<\nimfdFlap takes either \"prepare\", \"build\" or \"wcgw\", which is my version of help because I think it's funny, as an argument. Chances are I ran into your problem setting up precice and openfoam in which case I will explain it there. \"imfdFlap build\" is automatically invoked by allrun.sh.\n\n>>Instructions:<<\nRun \"" + executableName + " prepare + <name of directory to contain all files (defaults to: " + workDir + ")>\" to generate installation scripts and start the installation with \"allrun.sh\"." << endl;
    }
}
