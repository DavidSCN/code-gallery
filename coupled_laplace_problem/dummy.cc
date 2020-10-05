// To compile use:
// mpic++ -I$PRECICE_ROOT/src main.cpp -lprecice -o solverdummy

#include <iostream>
#include <sstream>

#include "precice/SolverInterface.hpp"

int
main()
{
  int commRank = 0;
  int commSize = 1;

  using namespace precice;
  using namespace precice::constants;

  const std::string configFileName("precice-config.xml");
  const std::string solverName("dummy-participant");
  const std::string meshName("dummy-mesh");
  const std::string dataWriteName("boundary-data");
  const std::string dataReadName("dummy");

  std::cout << "DUMMY: Running solver dummy with preCICE config file \""
            << configFileName << "\", participant name \"" << solverName
            << "\", and mesh name \"" << meshName << "\".\n";

  SolverInterface interface(solverName, configFileName, commRank, commSize);

  int meshID           = interface.getMeshID(meshName);
  int dimensions       = interface.getDimensions();
  int numberOfVertices = 10;

  const int readDataID  = interface.getDataID(dataReadName, meshID);
  const int writeDataID = interface.getDataID(dataWriteName, meshID);

  std::vector<double> readData(numberOfVertices);
  std::vector<double> writeData(numberOfVertices);
  std::vector<int>    vertexIDs(numberOfVertices);
  std::vector<double> vertices(numberOfVertices * dimensions);

  double deltaX = 2. / (numberOfVertices - 1);
  for (int i = 0; i < numberOfVertices; ++i)
    {
      for (int j = 0; j < dimensions; ++j)
        {
          const unsigned int index = dimensions * i + j;
          if (j == 0)
            vertices[index] = 1;
          else
            vertices[index] = 1 - deltaX * i;
        }

      readData[i]  = 0;
      writeData[i] = 2;
    }
  //  for(uint i=0;i<vertices.size();++i)
  //    std::cout<<vertices[i]<<std::endl;

  interface.setMeshVertices(meshID,
                            numberOfVertices,
                            vertices.data(),
                            vertexIDs.data());


  double dt = interface.initialize();

  while (interface.isCouplingOngoing())
    {
      if (interface.isActionRequired(actionWriteIterationCheckpoint()))
        {
          std::cout << "DUMMY: Writing iteration checkpoint\n";
          interface.markActionFulfilled(actionWriteIterationCheckpoint());
        }

      if (interface.isReadDataAvailable())
        {
          interface.readBlockScalarData(readDataID,
                                        numberOfVertices,
                                        vertexIDs.data(),
                                        readData.data());
        }

      for (int i = 0; i < numberOfVertices; i++)
        {
          std::cout << readData[i] << std::endl;
        }

      if (interface.isWriteDataRequired(dt))
        {
          interface.writeBlockScalarData(writeDataID,
                                         numberOfVertices,
                                         vertexIDs.data(),
                                         writeData.data());
        }

      dt = interface.advance(dt);

      if (interface.isActionRequired(actionReadIterationCheckpoint()))
        {
          std::cout << "DUMMY: Reading iteration checkpoint\n";
          interface.markActionFulfilled(actionReadIterationCheckpoint());
        }
      else
        {
          std::cout << "DUMMY: Advancing in time\n";
        }
    }

  interface.finalize();
  std::cout << "DUMMY: Closing C++ solver dummy...\n";

  return 0;
}
