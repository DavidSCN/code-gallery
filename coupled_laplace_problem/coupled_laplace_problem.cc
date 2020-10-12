// The included deal.II header files are exactly the same as in the step-4
// tutorial program
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
// In addition to the deal.II header files, we include the preCICE API in order
// to obtain access to preCICE specific functionality
#include <precice/SolverInterface.hpp>

#include <fstream>
#include <iostream>

using namespace dealii;

// @sect3{Configuration parameters}
//
// Here, we set up a simple hard-coded Struct containing the configuration of
// our simulation setup. The configuration includes the name of the
// configuration file, the name of the simulation participants and the names of
// the exchanged data. The same namings are used again in the precice
// configuration file, so that preCICE knows, which data should be exchanged.
// For real application cases, the configuration is likely better handled by a
// parameter file.
struct CouplingParamters
{
  const std::string config_file      = "precice-config.xml";
  const std::string participant_name = "laplace-solver";
  const std::string mesh_name        = "original-mesh";
  const std::string write_data_name  = "dummy";
  const std::string read_data_name   = "boundary-data";
};


// @sect4{The Adapter class}
//
// The Adapter class keeps all functionalities to couple the deal.II solver code
// to other solvers with preCICE, i.e., data structures are set up and all
// relevant information is passed to preCICE.

template <int dim, typename ParameterClass>
class Adapter
{
public:
  Adapter(const ParameterClass &   parameters,
          const types::boundary_id deal_boundary_interface_id);

  ~Adapter();

  void
  initialize(const DoFHandler<dim> &                    dof_handler,
             std::map<types::global_dof_index, double> &boundary_data);

  void
  advance(std::map<types::global_dof_index, double> &boundary_data,
          const double                               computed_timestep_length);

  // public precice solverinterface
  precice::SolverInterface precice;

  // Boundary ID of the deal.II triangulation, associated with the coupling
  // interface. The variable is defined in the constructor of this class and
  // intentionally public so that it can be used during the grid generation and
  // system assembly. The only thing, one needs to make sure is, that this ID is
  // unique for a particular triangulation.
  const unsigned int deal_boundary_interface_id;

private:
  // preCICE related initializations
  // These variables are specified in and read from a parameter file, which is
  // in this simple tutorial program the CouplingParameter struct already
  // introduced in the beginning.
  const std::string mesh_name;
  const std::string read_data_name;
  const std::string write_data_name;

  // These IDs are filled by preCICE during the initialization. We set a default
  // value of -1 in order to detect potential errors more easily.
  int mesh_id           = -1;
  int read_data_id      = -1;
  int write_data_id     = -1;
  int n_interface_nodes = -1;

  // Dof IndexSet, containing relevant coupling dof indices at the coupling
  // boundary
  IndexSet coupling_dofs;

  // Data containers which are passed to preCICE in an appropriate preCICE
  // specific format
  std::vector<int>    interface_nodes_ids;
  std::vector<double> read_data;
  std::vector<double> write_data;


  void
  format_precice_to_deal(
    std::map<types::global_dof_index, double> &boundary_data) const;
};



// In the constructor of the Adapter class, we set up the preCICE
// Solverinterface. Here, we need to tell preCICE our name as participant of the
// simulation and the name of the preCICE-configuration file. Both have already
// been specified in the CouplingParameter class above. Thus, we pass the class
// directly to the constructor and read out all relevant information. As a
// second parameter, we need to specify the boundary ID of our triangulation,
// which is associated with the coupling interface.
template <int dim, typename ParameterClass>
Adapter<dim, ParameterClass>::Adapter(
  const ParameterClass &   parameters,
  const types::boundary_id deal_boundary_interface_id)
  : precice(parameters.participant_name,
            parameters.config_file,
            Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),
            Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
  , deal_boundary_interface_id(deal_boundary_interface_id)
  , mesh_name(parameters.mesh_name)
  , read_data_name(parameters.read_data_name)
  , write_data_name(parameters.write_data_name)
{}



// The destructor of the Adapter class finalizes the Solverinterface, so that
// the communication is terminated and memory is deallocated.
template <int dim, typename ParameterClass>
Adapter<dim, ParameterClass>::~Adapter()
{
  precice.finalize();
}



// This function initializes preCICE (e.g. establishes communication channels
// and allocates memory) and passes all relevant data to preCICE. For surface
// coupling, relevant data is especially the location of the data points at the
// assoociated interface(s). The `boundary_data` is an empty map, which is
// filled by preCICE, i.e., the information of our dummy participant. Throughout
// the system assembly, the map can directly be used in order to apply the
// Dirichlet boundary conditions in the linear system.
template <int dim, typename ParameterClass>
void
Adapter<dim, ParameterClass>::initialize(
  const DoFHandler<dim> &                    dof_handler,
  std::map<types::global_dof_index, double> &boundary_data)
{
  Assert(dim > 1, ExcNotImplemented());
  Assert(dim == precice.getDimensions(),
         ExcDimensionMismatch(dim, precice.getDimensions()));

  // In a first step, we get precice specific IDs from precice and store them in
  // the respective variables. Later, they are used for data transfer.
  mesh_id = precice.getMeshID(mesh_name);

  if (precice.hasData(read_data_name, mesh_id))
    read_data_id = precice.getDataID(read_data_name, mesh_id);

  if (precice.hasData(write_data_name, mesh_id))
    write_data_id = precice.getDataID(write_data_name, mesh_id);

  // Afterwards, we extract the number of interface nodes and the coupling DoFs
  // at the coupling interface from our deal.II solver via
  // `extract_boundary_dofs()`
  std::set<types::boundary_id> couplingBoundary;
  couplingBoundary.insert(deal_boundary_interface_id);

  // The `ComponentMask()` might be important in case we deal with vector valued
  // problems, because vector valued problems have a DoF for each component.
  DoFTools::extract_boundary_dofs(dof_handler,
                                  ComponentMask(),
                                  coupling_dofs,
                                  couplingBoundary);

  // The coupling DoFs are used to set up the `boundary_data` map. At the end,
  // we associate here each DoF with a respective boundary value.
  for (auto i : coupling_dofs)
    boundary_data.insert(std::pair<types::global_dof_index, double>(i, 0));

  // Since we deal with a scalar problem, the number of DoFs at the particular
  // interface corresponds to the number of interface nodes.
  n_interface_nodes = coupling_dofs.n_elements();

  std::cout << "\t Number of coupling nodes:     " << n_interface_nodes
            << std::endl;

  // Now, we need to tell preCICE the coordinates of the interface nodes. The
  // preCICE API requires the utilization of std::vectors. Hence, we set up a
  // vector to pass the node positions to preCICE. Each node is specified only
  // once.
  std::vector<double> interface_nodes_positions(dim * n_interface_nodes);

  // Set up the appropriate size of the data container needed for data
  // exchange. Here, we deal with a scalar problem, so that only a scalar value
  // is read/written per interface node.
  write_data.resize(n_interface_nodes);
  read_data.resize(n_interface_nodes);
  // The IDs are again filled by preCICE during the initializations.
  interface_nodes_ids.resize(n_interface_nodes);

  // The node location is obtained using `map_dofs_to_support_points()`.
  std::map<types::global_dof_index, Point<dim>> support_points;

  // We use here a simple Q1 mapping. In case one has more complex
  // geomtries, you might want to change this to a higher order mapping.
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(),
                                       dof_handler,
                                       support_points);

  // `support_points` contains now the coordinates of all DoFs. In the next
  // step, the relevant coordinates are extracted using the IndexSet with the
  // extracted coupling_dofs.

  // preCICE expects all data in the format [x0, y0, z0, x1, y1 ...]. In 2D, the
  // last z coordinate needs to be omitted.
  int node_position_iterator = 0;
  for (auto element : coupling_dofs)
    {
      for (int i = 0; i < dim; ++i)
        interface_nodes_positions[node_position_iterator * dim + i] =
          support_points[element][i];

      ++node_position_iterator;
    }

  // Now we have all information to define the coupling mesh and pass the
  // information to preCICE.
  precice.setMeshVertices(mesh_id,
                          n_interface_nodes,
                          interface_nodes_positions.data(),
                          interface_nodes_ids.data());

  // Then, we initialize preCICE internally calling the API function
  // `initialize()`
  precice.initialize();

  // write initial writeData to preCICE if required
  // TODO: Should we remove this block here? No data is written, but it would
  // illustrate the concept.
  if (precice.isActionRequired(precice::constants::actionWriteInitialData()))
    {
      // store initial write_data for precice in write_data
      // Comment about extracting data from deal to precice
      precice.writeBlockScalarData(write_data_id,
                                   n_interface_nodes,
                                   interface_nodes_ids.data(),
                                   write_data.data());

      precice.markActionFulfilled(precice::constants::actionWriteInitialData());
    }

  // TODO: Discuss position of this function: inside the if statement leads to a
  // precice error, because no initial data is required, but it needs to be
  // initialized
  precice.initializeData();

  // read initial readData from preCICE if required for the first time step
  // FIXME: Data is already exchanged here
  if (precice.isReadDataAvailable())
    {
      precice.readBlockScalarData(read_data_id,
                                  n_interface_nodes,
                                  interface_nodes_ids.data(),
                                  read_data.data());

      // This is the opposite direction as above. See comment there.
      format_precice_to_deal(boundary_data);
    }
}


// The function `advance()` is called in the main time loop after the
// computation in each time step in the individual solver has finished. Here,
// coupling data is passed to preCICE and obtained from other participants. In
// case of the simplified unidirectional coupling, we just obtain data from our
// dummy participant.
template <int dim, typename ParameterClass>
void
Adapter<dim, ParameterClass>::advance(
  std::map<types::global_dof_index, double> &boundary_data,
  const double                               computed_timestep_length)
{
  // This is essentially the same as during initialization. We have already all
  // IDs and just need to convert our obtained data to the preCICE compatible
  // `write_data` vector. Since we deal here with a unidirectional example,
  // where the Laplace solver doesn't write any data, we don't need to transform
  // any data and the following code block is not accessed, i.e.
  // `precice.hasData()` will return false.
  // TODO: Should we remove this block here? No data is written, but it would
  // illustrate the concept.
  if (precice.hasData(write_data_name, mesh_id))
    if (precice.isWriteDataRequired(computed_timestep_length))
      {
        precice.writeBlockScalarData(write_data_id,
                                     n_interface_nodes,
                                     interface_nodes_ids.data(),
                                     write_data.data());
      }

  // Here, we specify the computed time step length and pass it to preCICE
  precice.advance(computed_timestep_length);

  // As a next step, we obtain data, i.e. the boundary condition, from another
  // participant. We have already all IDs and just need to convert our obtained
  // data to the deal.II compatible 'boundary map' , which is done in the
  // format_deal_to_precice function. All this is of course only done in
  // case write data is required.
  if (precice.isReadDataAvailable())
    {
      precice.readBlockScalarData(read_data_id,
                                  n_interface_nodes,
                                  interface_nodes_ids.data(),
                                  read_data.data());

      format_precice_to_deal(boundary_data);
    }
}



// This function takes the std::vector obtained by preCICE in `read_data` and
// inserts the values to the right position in the boundary map used throughout
// our deal.II solver for Dirichlet boundary conditions. The function is only
// used internally in the Adapter class and not called in the solver itself. The
// order, in which preCICE sorts the data in the `read_data` vector is exactly
// the same as the order of the initially passed vertices coordinates.
template <int dim, typename ParameterClass>
void
Adapter<dim, ParameterClass>::format_precice_to_deal(
  std::map<types::global_dof_index, double> &boundary_data) const
{
  // We stored already the coupling DoF indices in the `boundary_data` map, so
  // that we can simply iterate over all keys in the map.
  auto dof_component = boundary_data.begin();
  for (int i = 0; i < n_interface_nodes; ++i)
    {
      AssertIndexRange(i, read_data.size());
      boundary_data[dof_component->first] = read_data[i];
      ++dof_component;
    }
}


// The solver class is essentially the same as in step-4. Comments are added at
// any point, where the workflow is different due to the coupling.
template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem();

  void
  run();

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results() const;

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;

  // Here, we allocate all structures required for the preCICE coupling. The map
  // is used to apply Dirichlet boundary conditions and filled in the Adapter
  // class with data from our dummy participant. The CouplingParameters hold the
  // preCICE configuration as described above. The interdace boundary ID is the
  // ID associated to our coupling interface and needs to be specified, when we
  // set up the Adapter class object, because we pass it directly to the
  // Constructor of this class.
  std::map<types::global_dof_index, double> boundary_data;
  CouplingParamters                         parameters;
  const types::boundary_id                  interface_boundary_id;
  Adapter<dim, CouplingParamters>           adapter;
};



template <int dim>
class RightHandSide : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};

template <int dim>
double
RightHandSide<dim>::value(const Point<dim> &p,
                          const unsigned int /*component*/) const
{
  double return_value = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value += 4.0 * std::pow(p(i), 4.0);

  return return_value;
}


template <int dim>
double
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  return p.square();
}



template <int dim>
LaplaceProblem<dim>::LaplaceProblem()
  : fe(1)
  , dof_handler(triangulation)
  , interface_boundary_id(1)
  , adapter(parameters, interface_boundary_id)
{}


template <int dim>
void
LaplaceProblem<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(4);

  for (const auto &cell : triangulation.active_cell_iterators())
    for (auto f : GeometryInfo<dim>::face_indices())
      {
        const auto face = cell->face(f);

        // We choose (arbitrarily) the boundary in positive x direction for the
        // interface coupling.
        if (face->at_boundary() && (face->center()[0] == 1))
          face->set_boundary_id(interface_boundary_id);
      }

  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}


template <int dim>
void
LaplaceProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}



template <int dim>
void
LaplaceProblem<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx

            const auto x_q = fe_values.quadrature_point(q_index);
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            right_hand_side.value(x_q) *        // f(x_q)
                            fe_values.JxW(q_index));            // dx
          }

      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        {
          for (const unsigned int j : fe_values.dof_indices())
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }
  {
    // At first, we apply the Dirichlet boundary condition from step-4, as
    // usual.
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
  }
  {
    // Afterwards, we apply the coupling boundary condition. The `boundary_data`
    // has already been filled by preCICE.
    MatrixTools::apply_boundary_values(boundary_data,
                                       system_matrix,
                                       solution,
                                       system_rhs);
  }
}



template <int dim>
void
LaplaceProblem<dim>::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}



template <int dim>
void
LaplaceProblem<dim>::output_results() const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk");
  data_out.write_vtk(output);
}



template <int dim>
void
LaplaceProblem<dim>::run()
{
  std::cout << "Solving problem in " << dim << " space dimensions."
            << std::endl;

  make_grid();
  setup_system();

  // After we set up out system, we initialize preCICE using the functionalities
  // of our Adapter.
  adapter.initialize(dof_handler, boundary_data);

  // preCICE steers the coupled simulation completely. The steering is provided
  // by the function `isCouplingOngoing`, which tells the solver, whether the
  // time loop has already been finished or not. This tutorial solves 'just' an
  // explicit coupled stationary problem, so that the return of this statement
  // is trivial, but it is kept here to illustrate the general concept.
  while (adapter.precice.isCouplingOngoing())
    {
      // In the time loop, we assemble the coupled system and solve it, as
      // usual. According to our configuration, we obtained already the data of
      // our dummy participant during the initialization.
      assemble_system();
      solve();

      // After we solved the system, we advance with our system to the next time
      // level. In a coupled simulation, we would pass our calculated data to
      // preCICE and obtain data from other participants. Here, we simply obtain
      // data from the dummy participant. Since we wish to solve a stationary
      // problem, we set the time-step size equal to the time-window size of the
      // coupled system, which is 1.
      adapter.advance(boundary_data, 1);

      // In case our time step has been completed, we write the results to an
      // output file.
      if (adapter.precice.isTimeWindowComplete())
        output_results();
    }
}



int
main(int argc, char **argv)
{
  // TODO: Should we keep the MPI dependency here? It makes things more clear in
  // the initialization of the Adapter, but requires DEAL_II_WITH_MPI to be
  // turned on I guess.
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  LaplaceProblem<2> laplace_problem;
  laplace_problem.run();

  return 0;
}
