
#include <iostream>
#include <sstream>

#include "precice/SolverInterface.hpp"

int
main()
{
  const int commRank = 0;
  const int commSize = 1;

  using namespace precice;
  using namespace precice::constants;

  // Configuration
  const std::string configFileName("precice-config.xml");
  const std::string solverName("dummy-participant");
  const std::string meshName("dummy-mesh");
  const std::string dataWriteName("boundary-data");
  const std::string dataReadName("dummy");

  std::cout << "DUMMY: Running solver dummy with preCICE config file \""
            << configFileName << "\", participant name \"" << solverName
            << "\", and mesh name \"" << meshName << "\".\n";

  SolverInterface interface(solverName, configFileName, commRank, commSize);

  const int meshID           = interface.getMeshID(meshName);
  const int dimensions       = interface.getDimensions();
  const int numberOfVertices = 6;


  const int readDataID = interface.hasData(dataReadName, meshID) ?
                           interface.getDataID(dataReadName, meshID) :
                           -1;
  const int writeDataID = interface.hasData(dataWriteName, meshID) ?
                            interface.getDataID(dataWriteName, meshID) :
                            -1;

  // Set up data structures
  std::vector<double> readData(numberOfVertices);
  std::vector<double> writeData(numberOfVertices);
  std::vector<int>    vertexIDs(numberOfVertices);
  std::vector<double> vertices(numberOfVertices * dimensions);

  // Define a boundary mesh
  const double deltaX = 2. / (numberOfVertices - 1);
  for (int i = 0; i < numberOfVertices; ++i)
    {
      for (int j = 0; j < dimensions; ++j)
        {
          const unsigned int index = dimensions * i + j;
          // The x-coordinate is always 1, i.e., the boundary is parallel to the
          // y-axis. The y-coordinate is descending from 1 to -1.
          if (j == 0)
            vertices[index] = 1;
          else
            vertices[index] = 1 - deltaX * i;
        }

      // Specify the actual data we want to pass to the other participant. We
      // enforce continuity to adjacent boundaries, where the value is set to 2
      // by specifying the first and the last value. The remaining vertices are
      // associated with a random value between 0 and 2.
      if (i == 0)
        writeData[i] = 2;
      else if (i == numberOfVertices - 1)
        writeData[i] = 2;
      else
        writeData[i] = (double)(rand() % 20) / 10.;
    }

  // Pass the vertices to preCICE
  interface.setMeshVertices(meshID,
                            numberOfVertices,
                            vertices.data(),
                            vertexIDs.data());

  // initialize the Solverinterface
  double dt = interface.initialize();

  if (interface.isActionRequired(actionWriteInitialData()))
    {
      std::cout << "DUMMY: Writing initial data \n";
      interface.writeBlockScalarData(writeDataID,
                                     numberOfVertices,
                                     vertexIDs.data(),
                                     writeData.data());

      interface.markActionFulfilled(actionWriteInitialData());

      interface.initializeData();
    }

  while (interface.isCouplingOngoing())
    {
      if (interface.isWriteDataRequired(dt))
        {
          std::cout << "DUMMY: Writing coupling data \n";
          interface.writeBlockScalarData(writeDataID,
                                         numberOfVertices,
                                         vertexIDs.data(),
                                         writeData.data());
        }

      dt = interface.advance(dt);
      std::cout << "DUMMY: Advancing in time\n";

      // FIXME: Should retun false
      if (interface.isReadDataAvailable() && false)
        {
          std::cout << "DUMMY: Reading coupling data \n";
          interface.readBlockScalarData(readDataID,
                                        numberOfVertices,
                                        vertexIDs.data(),
                                        readData.data());
        }
    }


  interface.finalize();
  std::cout << "DUMMY: Closing C++ solver dummy...\n";

  return 0;
}
