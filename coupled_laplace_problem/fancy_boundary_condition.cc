
#include <iostream>
#include <sstream>

#include "precice/SolverInterface.hpp"

// Function to generate boundary values in each time step
void
define_boundary_values(std::vector<double> &boundary_data,
                       const double         time,
                       const double         end_time)
{
  // Scale the current time value
  const double relative_time = time / end_time;
  // Define the amplitude. Values run from -0.5 to 0.5
  const double amplitude = (relative_time - 0.5);
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
  std::cout << "Boundary participant: starting... \n";

  // Adjust to MPI rank and size for parallel computation
  const int commRank = 0;
  const int commSize = 1;

  using namespace precice;

  // Configuration
  const std::string configFileName("precice-config.xml");
  const std::string solverName("boundary-participant");
  const std::string meshName("boundary-mesh");
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
  std::cout << "Boundary participant: defining boundary mesh \n";
  const double length = 2;
  const double xCoord = 1;
  const double deltaY = length / (numberOfVertices - 1);
  for (int i = 0; i < numberOfVertices; ++i)
    for (int j = 0; j < dimensions; ++j)
      {
        const unsigned int index = dimensions * i + j;
        // The x-coordinate is always 1, i.e., the boundary is parallel to the
        // y-axis. The y-coordinate is descending from 1 to -1.
        if (j == 0)
          vertices[index] = xCoord;
        else
          vertices[index] = 1 - deltaY * i;
      }

  // Pass the vertices to preCICE
  precice.setMeshVertices(meshID,
                          numberOfVertices,
                          vertices.data(),
                          vertexIDs.data());

  // initialize the Solverinterface
  double dt = precice.initialize();

  // Start time loop
  const double end_time = 1;
  double       time     = 0;
  while (precice.isCouplingOngoing())
    {
      // Generate new boundary data
      define_boundary_values(writeData, time, end_time);

      {
        std::cout << "Boundary participant: writing coupling data \n";
        precice.writeBlockScalarData(dataID,
                                     numberOfVertices,
                                     vertexIDs.data(),
                                     writeData.data());
      }

      dt = precice.advance(dt);
      std::cout << "Boundary participant: advancing in time\n";

      time += dt;
    }

  std::cout << "Boundary participant: closing...\n";

  return 0;
}
