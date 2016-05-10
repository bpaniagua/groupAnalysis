/*
* DiagnosticIndex research code
* Beatriz Paniagua UNC
*
*/
#include "PCAModelBuilder.h"
#include "StatisticalModel.h"
#include "DataManager.h"
#include "vtkStandardMeshRepresenter.h"
#include <boost/scoped_ptr.hpp>

#include "vtkDirectory.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "CommonTypes.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include <vtkDelimitedTextReader.h> //CSVloader
#include <vtkTable.h> //CSVloader

using namespace statismo;


typedef vtkStandardMeshRepresenter RepresenterType;
typedef DataManager<vtkPolyData> DataManagerType;
typedef PCAModelBuilder<vtkPolyData> ModelBuilderType;
typedef StatisticalModel<vtkPolyData> StatisticalModelType;
typedef std::vector<std::string> StringVectorType;

int getdir (std::string dir, std::vector<std::string> &files, const std::string& extension=".*") {
vtkDirectory *directory = vtkDirectory::New();
directory->Open(dir.c_str());
for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
const char* filename = directory->GetFile(i);
if (extension == ".*" || std::string(filename).find(extension) != std::string::npos)
files.push_back(filename);
}
directory->Delete();
return 0;
}

vtkPolyData* loadVTKPolyData(const std::string& filename) {
vtkPolyDataReader* reader = vtkPolyDataReader::New();
reader->SetFileName(filename.c_str());
reader->Update();
vtkPolyData* pd = vtkPolyData::New();
pd->ShallowCopy(reader->GetOutput());
return pd;
}

void saveSample(const vtkPolyData* pd, const std::string& resdir, const std::string& basename) {
    std::string filename = resdir +std::string("/") + basename;

    vtkPolyDataWriter* w = vtkPolyDataWriter::New();
#if (VTK_MAJOR_VERSION == 5 )
    w->SetInput(const_cast<vtkPolyData*>(pd));
#else
    w->SetInputData(const_cast<vtkPolyData*>(pd));
#endif
    w->SetFileName(filename.c_str());
    w->Update();
}

vtkPoints* subtractMesh (vtkPolyData* mesh1, vtkPolyData* mesh2)
{
    vtkPoints* subtractedPoints = vtkPoints::New();

    vtkPoints* pointsMesh1 = vtkPoints::New();
    pointsMesh1 = mesh1->GetPoints();
    vtkPoints* pointsMesh2 = vtkPoints::New();
    pointsMesh2 = mesh2->GetPoints();

    for( unsigned int pointID = 0; pointID < mesh1->GetNumberOfPoints(); pointID++ )
    {
        double mesh1point[3];
        double mesh2point[3];
        double tmpPoint[3];

        pointsMesh1->GetPoint(pointID, mesh1point);
        pointsMesh2->GetPoint(pointID, mesh2point);

        for( unsigned int dim = 0; dim < 3; dim++ )
        {
            tmpPoint[dim] = mesh1point[dim] - mesh2point[dim];
        }
        subtractedPoints->InsertPoint(pointID, tmpPoint);
    }

    return subtractedPoints;
}

bool is_file_exist(std::string fileName)
{
    std::ifstream infile(fileName.c_str());
    return infile.good();
}

StatisticalModelType* buildSSM (std::string datadir, std::string filenamereference)
{
    StringVectorType filenames;
    getdir(datadir, filenames, ".vtk");
    if (filenames.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadir << " exiting.";
    exit(-1);
    }

    std::cout << "Building SSM of data in " << datadir << std::endl;
    vtkPolyData* reference = loadVTKPolyData(filenamereference);
    boost::scoped_ptr<RepresenterType> representer(RepresenterType::Create(reference));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManager(DataManagerType::Create(representer.get()));

    for (unsigned i = 0; i < filenames.size() ; i++)
    {
        vtkPolyData* dataset = loadVTKPolyData(datadir + "/" + filenames[i]);
        // We provde the filename as a second argument.
        // It will be written as metadata, and allows us to more easily figure out what we did later.
        dataManager->AddDataset(dataset, filenames[i]);
        // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
        dataset->Delete();
    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    ModelBuilderType* modelBuilder = ModelBuilderType::Create();
    StatisticalModelType* model = modelBuilder->BuildNewModel(dataManager->GetData(), 0.01);

    std::cout << "Total variance " << model->GetPCAVarianceVector().sum() << " number of Eigenmodes " << model->GetPCAVarianceVector().size() << std::endl;
    return model;
}

StatisticalModelType* buildSSMCrossvalidation (std::string datadir, std::string filenamereference, std::string filename)
{
    StringVectorType filenames;
    getdir(datadir, filenames, ".vtk");
    if (filenames.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadir << " exiting.";
    exit(-1);
    }

  //  std::cout << "Building SSM of data in " << datadir << std::endl;
    vtkPolyData* reference = loadVTKPolyData(filenamereference);
    boost::scoped_ptr<RepresenterType> representer(RepresenterType::Create(reference));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManager(DataManagerType::Create(representer.get()));

    for (unsigned i = 0; i < filenames.size() ; i++)
    {
        if (filename.compare(filenames[i]) != 0)
        {
            vtkPolyData* dataset = loadVTKPolyData(datadir + "/" + filenames[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManager->AddDataset(dataset, filenames[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
        }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    ModelBuilderType* modelBuilder = ModelBuilderType::Create();
    StatisticalModelType* model = modelBuilder->BuildNewModel(dataManager->GetData(), 0.01);

  //  std::cout << "Total variance " << model->GetPCAVarianceVector().sum() << " number of Eigenmodes " << model->GetPCAVarianceVector().size() << std::endl;
    return model;
}

void analyzeVariance (StatisticalModelType* model)
{
    VectorType PCAVariance = model->GetPCAVarianceVector();
    double totalVar = model->GetPCAVarianceVector().sum();
    std::cout << "### Model variance analysis , total variance " << totalVar << "###"<< std::endl;

    int sum =0;
    double per = 0;
    for (int i = 0 ; i < PCAVariance.size() ; i++)
    {
        sum = sum + PCAVariance[i];
        double percentage = (PCAVariance[i] * 100) / totalVar;
        per = per + percentage;
        std::cout << "Eigenvalue " << i << " variance " << PCAVariance[i] << " (" << percentage <<" / " << per << ")" << std::endl;
    }
}

double ComputeOAindex (vtkPoints* diff)
{

    double OAindex;
    double sum = 0;

    for( unsigned int pointID = 0; pointID < diff->GetNumberOfPoints(); pointID++ )
    {
        double point[3];

        diff->GetPoint(pointID, point);

        for( unsigned int dim = 0; dim < 3; dim++ )
        {
            sum = pow ( point[dim],2) + sum ;
        }

    }

    OAindex = sum / (diff->GetNumberOfPoints()*3);
    return OAindex;
}

int main (int argc, char ** argv)
{

    // Variables that can be changed by the user This will be implemented in the GUI
    int NumOfGroups = 8;

    std::string GroupLUT("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/LookUpGroup2.csv");
    std::string datadirModels("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/HD5Groups");

    std::string datadirHC("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/ControlGroup");
    std::string datadirOA("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/OA");
    std::string datadirBoth("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Both");

    std::string datadirG01("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/G01");
    std::string datadirG02("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/G02");
    std::string datadirG03("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/G03");
    std::string datadirG04("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/G04");
    std::string datadirG05("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/G05");
    std::string datadirG06("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/G06");
    std::string datadirG07("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/G07");

    std::string datadirRepresenters("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Representers");

    //Read .CSV look up table with group information and VTK
    vtkSmartPointer<vtkDelimitedTextReader> CSVreader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    CSVreader->SetFieldDelimiterCharacters(",");
    CSVreader->SetFileName(GroupLUT.c_str());
    CSVreader->SetHaveHeaders(false);
    CSVreader->Update();
    vtkSmartPointer<vtkTable> table = CSVreader->GetOutput();


    StringVectorType filenamesToClassify;
    getdir(datadirBoth, filenamesToClassify, ".vtk");
    if (filenamesToClassify.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirBoth << " exiting.";
    exit(-1);
    }

    std::cout << "Number of samples to classify " << filenamesToClassify.size() << std::endl;

    VectorType totalGroupVariances (8);
    double totalVarPooled=0;

    Eigen::MatrixXd OAindex_all (NumOfGroups,filenamesToClassify.size());
    VectorType OAindex (filenamesToClassify.size());
    VectorType OAindex_groupAssignment (filenamesToClassify.size());
    VectorType OAindex_groupReal (filenamesToClassify.size());

    for (unsigned sample = 0 ; sample < filenamesToClassify.size() ; sample++)
    {

        int group = 0;

        for ( int r = 0; r < table->GetNumberOfRows() ; r++ )
        {

            if ( filenamesToClassify[sample] == table->GetValue(r,0).ToString())
            {
                group = atoi(table->GetValue(r,1).ToString().c_str());
            }

         }

        std::cout << sample <<  ": Group assignment for " << filenamesToClassify[sample] << " is " << group << std::endl;
        OAindex_groupReal[sample] = group;


        //Computing variances per group and pooled
    //    totalGroupVariances(0) = modelG01->GetPCAVarianceVector().sum();
    //    totalGroupVariances(1) = modelG02->GetPCAVarianceVector().sum();
     //   totalGroupVariances(2) = modelG03->GetPCAVarianceVector().sum();
     //   totalGroupVariances(3) = modelG04->GetPCAVarianceVector().sum();
     //   totalGroupVariances(4) = modelG05->GetPCAVarianceVector().sum();
      //  totalGroupVariances(5) = modelG06->GetPCAVarianceVector().sum();
       // totalGroupVariances(6) = modelG07->GetPCAVarianceVector().sum();
      //  totalGroupVariances(7) = modelHC->GetPCAVarianceVector().sum();

        double sum = 0;
        for (int i = 0; i <totalGroupVariances.size(); i++)
        {
            sum = sum + totalGroupVariances(i);
        }
        totalVarPooled = sum/totalGroupVariances.size();



        // ###################################
        //Computing OAindex for all samples
        // ###################################

        // Read ShapeOA in vtk format
        std::string ShapeOAFilename = datadirBoth + "/" + filenamesToClassify[sample];
        vtkPolyData* VTKShapeOA = loadVTKPolyData(ShapeOAFilename);

        // Load all avgs from each group
        std::string filenameAvgG01 = datadirRepresenters + "/avgG01.vtk";
        vtkPolyData* VTKavgG01 = loadVTKPolyData(filenameAvgG01);
        std::string filenameAvgG02 = datadirRepresenters + "/avgG02.vtk";
        vtkPolyData* VTKavgG02 = loadVTKPolyData(filenameAvgG02);
        std::string filenameAvgG03 = datadirRepresenters + "/avgG03.vtk";
        vtkPolyData* VTKavgG03 = loadVTKPolyData(filenameAvgG03);
        std::string filenameAvgG04 = datadirRepresenters + "/avgG04.vtk";
        vtkPolyData* VTKavgG04 = loadVTKPolyData(filenameAvgG04);
        std::string filenameAvgG05 = datadirRepresenters + "/avgG05.vtk";
        vtkPolyData* VTKavgG05 = loadVTKPolyData(filenameAvgG05);
        std::string filenameAvgG06 = datadirRepresenters + "/avgG06.vtk";
        vtkPolyData* VTKavgG06 = loadVTKPolyData(filenameAvgG06);
        std::string filenameAvgG07 = datadirRepresenters + "/avgG07.vtk";
        vtkPolyData* VTKavgG07 = loadVTKPolyData(filenameAvgG07);
        std::string filenameAvgHC = datadirRepresenters + "/avgHC.vtk";
        vtkPolyData* VTKavgHC = loadVTKPolyData(filenameAvgHC);

        vtkPoints* diffG01 = subtractMesh(VTKShapeOA,VTKavgG01);
        vtkPoints* diffG02 = subtractMesh(VTKShapeOA,VTKavgG02);
        vtkPoints* diffG03 = subtractMesh(VTKShapeOA,VTKavgG03);
        vtkPoints* diffG04 = subtractMesh(VTKShapeOA,VTKavgG04);
        vtkPoints* diffG05 = subtractMesh(VTKShapeOA,VTKavgG05);
        vtkPoints* diffG06 = subtractMesh(VTKShapeOA,VTKavgG06);
        vtkPoints* diffG07 = subtractMesh(VTKShapeOA,VTKavgG07);
        vtkPoints* diffHC = subtractMesh(VTKShapeOA,VTKavgHC);

        double minOAindex = 99999999999999;
        double maxOAindex = -1;
        double groupAssignment = 0;

        OAindex_all(0,sample) = ComputeOAindex(diffG01);
        OAindex_all(1,sample) = ComputeOAindex(diffG02);
        OAindex_all(2,sample) = ComputeOAindex(diffG03);
        OAindex_all(3,sample) = ComputeOAindex(diffG04);
        OAindex_all(4,sample) = ComputeOAindex(diffG05);
        OAindex_all(5,sample) = ComputeOAindex(diffG06);
        OAindex_all(6,sample) = ComputeOAindex(diffG07);
        OAindex_all(7,sample) = ComputeOAindex(diffHC);


     //   OAindex_all(0,sample) = modelG01->ComputeProbabilityOfDataset(VTKShapeOA);
      //  OAindex_all(1,sample) = modelG02->ComputeProbabilityOfDataset(VTKShapeOA);
      //  OAindex_all(2,sample) = modelG03->ComputeProbabilityOfDataset(VTKShapeOA);
       // OAindex_all(3,sample) = modelG04->ComputeProbabilityOfDataset(VTKShapeOA);
       // OAindex_all(4,sample) = modelG05->ComputeProbabilityOfDataset(VTKShapeOA);
      //  OAindex_all(5,sample) = modelG06->ComputeProbabilityOfDataset(VTKShapeOA);
     //   OAindex_all(6,sample) = modelG07->ComputeProbabilityOfDataset(VTKShapeOA);
     //   OAindex_all(7,sample) = modelHC->ComputeProbabilityOfDataset(VTKShapeOA);

        for (int i = 0 ; i < NumOfGroups ; i ++ )
        {
            if (OAindex_all(i,sample) < minOAindex)
            {
                minOAindex = OAindex_all(i,sample);
                groupAssignment = i+1;
            }

        }

        OAindex(sample) = minOAindex;
        OAindex_groupAssignment(sample) = groupAssignment;
        std::cout << "--------------------" << std::endl;

    }

    std::cout << "----- Writing OA index results -----" << std::endl;
//
    std::ofstream outputOAindexAll;
    outputOAindexAll.open( "/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/avgAnalysis/OAindexAll_groups.csv" );
    outputOAindexAll << "GroupID,";
    for (int ij = 0; ij < filenamesToClassify.size()-1; ij ++)
        outputOAindexAll << filenamesToClassify[ij] << "," ;
    outputOAindexAll << filenamesToClassify[filenamesToClassify.size()-1] << std::endl;

    for (int ij = 0; ij < OAindex_all.rows(); ij ++)
    {
        outputOAindexAll << ij+1 << "," ;
        for (int kl = 0; kl < OAindex_all.cols()-1; kl ++)
        {
            outputOAindexAll << OAindex_all(ij,kl) << "," ;
        }
        outputOAindexAll << OAindex_all(ij,(OAindex_all.cols()-1)) << std::endl;
    }
    outputOAindexAll.close();

    std::ofstream outputOAindex;
    outputOAindex.open( "/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/avgAnalysis/OAindex_groups.csv" );
    outputOAindex << "subjectID,OAindex,GroupAssignment,GroupReal,Classification" << std::endl;

    int missclassifications = 0;
    for (int ij = 0; ij < OAindex.size(); ij ++)
    {
        int missclassified = 0;
        if (OAindex_groupAssignment[ij] != OAindex_groupReal[ij])
        {
            missclassified = 1;
            missclassifications++;
        }
        outputOAindex << filenamesToClassify[ij] << "," << OAindex[ij] << "," << OAindex_groupAssignment[ij] << "," << OAindex_groupReal[ij] << "," << missclassified << std::endl;
    }
    outputOAindex << "Number of missclassifications " << missclassifications << "/" << (missclassifications * 100) / filenamesToClassify.size() << "%" << std::endl;
    outputOAindex.close();


}
