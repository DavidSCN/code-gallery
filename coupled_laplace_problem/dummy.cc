
#include <iostream>
#include <sstream>

#include "precice/SolverInterface.hpp"

void
define_boundary_values(std::vector<double> &boundary_data, double amplitude)
{
  // Specify the actual data we want to pass to the other participant. Here, we
  // choose a parabola with boundary values 2 in order to enforce continuity
  // to adjacent boundaries.
  const unsigned int n_elements = boundary_data.size();
  for (uint i = 0; i < n_elements; ++i)
    boundary_data[i] = amplitude * (i * ((n_elements - 1) - i)) + 2;
}


int
main()
{
  // Adjust to MPI rank and size for parallel computation
  const int commRank = 0;
  const int commSize = 1;

  using namespace precice;

  // Configuration
  const std::string configFileName("precice-config.xml");
  const std::string solverName("dummy-participant");
  const std::string meshName("dummy-mesh");
  const std::string dataWriteName("boundary-data");

  SolverInterface precice(solverName, configFileName, commRank, commSize);

  const int meshID           = precice.getMeshID(meshName);
  const int dimensions       = precice.getDimensions();
  const int numberOfVertices = 6;


  const int dataID = precice.getDataID(dataWriteName, meshID);

  // Set up data structures
  std::vector<double> writeData(numberOfVertices);
  std::vector<int>    vertexIDs(numberOfVertices);
  std::vector<double> vertices(numberOfVertices * dimensions);

  // Define a boundary mesh
  const double length  = 2;
  const double x_coord = 1;
  const double deltaX  = length / (numberOfVertices - 1);
  for (int i = 0; i < numberOfVertices; ++i)
    for (int j = 0; j < dimensions; ++j)
      {
        const unsigned int index = dimensions * i + j;
        // The x-coordinate is always 1, i.e., the boundary is parallel to the
        // y-axis. The y-coordinate is descending from 1 to -1.
        if (j == 0)
          vertices[index] = x_coord;
        else
          vertices[index] = 1 - deltaX * i;
      }

  // Start with an amplitude of -1
  define_boundary_values(writeData, -1);

  // Pass the vertices to preCICE
  precice.setMeshVertices(meshID,
                          numberOfVertices,
                          vertices.data(),
                          vertexIDs.data());

  // initialize the Solverinterface
  double dt = precice.initialize();

  if (precice.isActionRequired(constants::actionWriteInitialData()))
    {
      std::cout << "DUMMY: Writing initial data \n";
      precice.writeBlockScalarData(dataID,
                                   numberOfVertices,
                                   vertexIDs.data(),
                                   writeData.data());

      precice.markActionFulfilled(constants::actionWriteInitialData());

      precice.initializeData();
    }

  double end_time = 10;
  double time     = -end_time / 2;
  while (precice.isCouplingOngoing())
    {
      time += dt;
      define_boundary_values(writeData, time / (end_time / 2));

      if (precice.isWriteDataRequired(dt))
        {
          std::cout << "DUMMY: Writing coupling data \n";
          precice.writeBlockScalarData(dataID,
                                       numberOfVertices,
                                       vertexIDs.data(),
                                       writeData.data());
        }

      dt = precice.advance(dt);
      std::cout << "DUMMY: Advancing in time\n";
    }

  std::cout << "DUMMY: Closing C++ solver dummy...\n";

  return 0;
}
