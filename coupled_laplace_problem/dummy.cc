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

  std::string configFileName("precice-config.xml");
  std::string solverName("dummy-participant");
  std::string meshName("dummy-mesh");

  std::cout << "DUMMY: Running solver dummy with preCICE config file \""
            << configFileName << "\", participant name \"" << solverName
            << "\", and mesh name \"" << meshName << "\".\n";

  SolverInterface interface(solverName, configFileName, commRank, commSize);

  int         meshID     = interface.getMeshID(meshName);
  int         dimensions = interface.getDimensions();
  std::string dataWriteName("dummy");
  std::string dataReadName("solution");
  int         numberOfVertices = 10;

  const int readDataID  = interface.getDataID(dataReadName, meshID);
  const int writeDataID = interface.getDataID(dataWriteName, meshID);

  std::vector<double> readData(numberOfVertices);
  std::vector<double> writeData(numberOfVertices);
  std::vector<double> vertices(numberOfVertices * dimensions);
  std::vector<int>    vertexIDs(numberOfVertices);

  double deltaX = 2. / (numberOfVertices - 1);
  for (int i = 0; i < numberOfVertices; i++)
    {
      for (int j = 0; j < dimensions; j++)
        {
          const unsigned int index = numberOfVertices * i + j;
          if ((index) % 2 == 0)
            vertices[index] = 1;
          else
            vertices[index] = 1 - deltaX * i;
        }

      readData[i]  = 0;
      writeData[i] = 0;
    }

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
          interface.readBlockVectorData(readDataID,
                                        numberOfVertices,
                                        vertexIDs.data(),
                                        readData.data());
        }

      for (int i = 0; i < numberOfVertices * dimensions; i++)
        {
          writeData[i] = readData[i] + 1;
        }

      if (interface.isWriteDataRequired(dt))
        {
          interface.writeBlockVectorData(writeDataID,
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
