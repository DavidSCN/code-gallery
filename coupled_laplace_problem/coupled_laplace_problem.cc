/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author:
 */

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

#include <precice/SolverInterface.hpp>

#include <fstream>
#include <iostream>

using namespace dealii;

struct CouplingParamters
{
  std::string config_file      = "precice-config.xml";
  std::string participant_name = "laplace-solver";
  std::string mesh_name        = "original-mesh";
  std::string write_data_name  = "solution";
  std::string read_data_name   = "dummy";
};



/**
 * The Adapter class keeps all functionalities to couple deal.II to other
 * solvers with preCICE i.e. data structures are set up, necessary information
 * is passed to preCICE etc.
 */
template <int dim, typename VectorType, typename ParameterClass>
class Adapter
{
public:
  /**
   * @brief      Constructor, which sets up the precice Solverinterface
   *
   * @param[in]  parameters Parameter class, which hold the data specified
   *             in the parameters.prm file
   * @param[in]  deal_boundary_interface_id Boundary ID of the triangulation,
   *             which is associated with the coupling interface.
   */
  Adapter(const ParameterClass &parameters,
          const unsigned int    deal_boundary_interface_id);

  /**
   * @brief      Destructor, which additionally finalizes preCICE
   *
   */
  ~Adapter();

  /**
   * @brief      Initializes preCICE and passes all relevant data to preCICE
   *
   * @param[in]  dof_handler Initialized dof_handler
   * @param[in]  deal_to_precice Data, which should be given to preCICE and
   *             exchanged with other participants. Wether this data is
   *             required already in the beginning depends on your
   *             individual configuration and preCICE determines it
   *             automatically. In many cases, this data will just represent
   *             your initial condition.
   * @param[out] precice_to_deal Data, which is received from preCICE/ from
   *             other participants. Wether this data is useful already in
   *             the beginning depends on your individual configuration and
   *             preCICE determines it automatically. In many cases, this
   *             data will just represent the initial condition of other
   *             participants.
   */
  void
  initialize(const DoFHandler<dim> &dof_handler,
             const VectorType &     deal_to_precice,
             VectorType &           precice_to_deal);

  /**
   * @brief      Advances preCICE after every timestep, converts data formats
   *             between preCICE and dealii
   *
   * @param[in]  deal_to_precice Same data as in @p initialize_precice() i.e.
   *             data, which should be given to preCICE after each time step
   *             and exchanged with other participants.
   * @param[out] precice_to_deal Same data as in @p initialize_precice() i.e.
   *             data, which is received from preCICE/other participants
   *             after each time step and exchanged with other participants.
   * @param[in]  computed_timestep_length Length of the timestep used by
   *             the solver.
   */
  void
  advance(const VectorType &deal_to_precice,
          VectorType &      precice_to_deal,
          const double      computed_timestep_length);

  /**
   * @brief public precice solverinterface
   */

  precice::SolverInterface precice;

  // Boundary ID of the deal.II mesh, associated with the coupling
  // interface. The variable is public and should be used during grid
  // generation, but is also involved during system assembly. The only thing,
  // one needs to make sure is, that this ID is not given to another part of
  // the boundary e.g. clamped one.
  const unsigned int deal_boundary_interface_id;

private:
  // preCICE related initializations
  // These variables are specified and read from the parameter file
  const std::string mesh_name;
  const std::string read_data_name;
  const std::string write_data_name;

  // These IDs are given by preCICE during initialization
  int mesh_id           = -1;
  int read_data_id      = -1;
  int write_data_id     = -1;
  int n_interface_nodes = -1;

  // Dof IndexSets of the global deal.II vectors, containing relevant
  // coupling dof indices
  IndexSet coupling_dofs;

  // Data containers which are passed to preCICE in an appropriate preCICE
  // specific format
  std::vector<int>    interface_nodes_ids;
  std::vector<double> read_data;
  std::vector<double> write_data;


  /**
   * @brief format_deal_to_precice Formats a global deal.II vector of type
   *        VectorType to a std::vector for preCICE. This functions is only
   *        used internally in the class and should not be called in the
   *        solver.
   *
   * @param[in] deal_to_precice Global deal.II vector of VectorType. The
   *            result (preCICE specific vector) is stored in the class in
   *            the variable 'write_data'.
   *
   * @note  The order, in which preCICE obtains data from the solver, needs
   *        to be consistent with the order of the initially passed vertices
   *        coordinates.
   */
  void
  format_deal_to_precice(const VectorType &deal_to_precice);

  /**
   * @brief format_precice_to_deal Takes the std::vector obtained by preCICE
   *        in 'read_data' and inserts the values to the right position in
   *        the global deal.II vector of size n_global_dofs. This is the
   *        opposite functionality as @p foramt_precice_to_deal(). This
   *        functions is only used internally in the class and should not
   *        be called in the solver.
   *
   * @param[out] precice_to_deal Global deal.II vector of VectorType and
   *             size n_global_dofs.
   *
   * @note  The order, in which preCICE obtains data from the solver, needs
   *        to be consistent with the order of the initially passed vertices
   *        coordinates.
   */
  void
  format_precice_to_deal(VectorType &precice_to_deal) const;
};



template <int dim, typename VectorType, typename ParameterClass>
Adapter<dim, VectorType, ParameterClass>::Adapter(
  const ParameterClass &parameters,
  const unsigned int    deal_boundary_interface_id)
  : precice(parameters.participant_name,
            parameters.config_file,
            Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),
            Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
  , deal_boundary_interface_id(deal_boundary_interface_id)
  , mesh_name(parameters.mesh_name)
  , read_data_name(parameters.read_data_name)
  , write_data_name(parameters.write_data_name)
{}



template <int dim, typename VectorType, typename ParameterClass>
Adapter<dim, VectorType, ParameterClass>::~Adapter()
{
  precice.finalize();
}



template <int dim, typename VectorType, typename ParameterClass>
void
Adapter<dim, VectorType, ParameterClass>::initialize(
  const DoFHandler<dim> &dof_handler,
  const VectorType &     deal_to_precice,
  VectorType &           precice_to_deal)
{
  Assert(dim > 1, ExcNotImplemented());
  Assert(dim == precice.getDimensions(), ExcInternalError());

  // get precice specific IDs from precice and store them in the class
  // they are later needed for data transfer
  mesh_id = precice.getMeshID(mesh_name);

  if (precice.hasData(read_data_name, mesh_id))
    read_data_id = precice.getDataID(read_data_name, mesh_id);

  if (precice.hasData(write_data_name, mesh_id))
    write_data_id = precice.getDataID(write_data_name, mesh_id);

  // get the number of interface nodes from deal.II
  // Therefore, we extract one component of the vector valued dofs and store
  // them in an IndexSet
  std::set<types::boundary_id> couplingBoundary;
  couplingBoundary.insert(deal_boundary_interface_id);

  DoFTools::extract_boundary_dofs(dof_handler,
                                  ComponentMask(),
                                  coupling_dofs,
                                  couplingBoundary);

  // Comment about scalar problems

  n_interface_nodes = coupling_dofs.n_elements();

  std::cout << "\t Number of coupling nodes:     " << n_interface_nodes
            << std::endl;

  // Set up a vector to pass the node positions to preCICE. Each node is
  // specified once. One needs to specify in the precice-config.xml, whether
  // the data is vector valued or not.
  std::vector<double> interface_nodes_positions(dim * n_interface_nodes);

  // Set up the appropriate size of the data container needed for data
  // exchange. Here, we deal with a vector valued problem for read and write
  // data namely displacement and forces.
  write_data.resize(n_interface_nodes, 0);
  read_data.resize(n_interface_nodes, 0);
  interface_nodes_ids.resize(n_interface_nodes);

  // get the coordinates of the interface nodes from deal.ii
  std::map<types::global_dof_index, Point<dim>> support_points;

  // We use here a simple Q1 mapping. In case one has more complex
  // geomtries, you might want to change this to a higher order mapping.
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(),
                                       dof_handler,
                                       support_points);

  // support_points contains now the coordinates of all dofs
  // in the next step, the relevant coordinates are extracted using the
  // IndexSet with the extracted coupling_dofs.

  // preCICE expects all data in the format [x0, y0, z0, x1, y1 ...]
  int node_position_iterator = 0;
  for (auto element : coupling_dofs)
    {
      for (int i = 0; i < dim; ++i)
        interface_nodes_positions[node_position_iterator * dim + i] =
          support_points[element][i];

      ++node_position_iterator;
    }

  // pass node coordinates to precice
  precice.setMeshVertices(mesh_id,
                          n_interface_nodes,
                          interface_nodes_positions.data(),
                          interface_nodes_ids.data());

  // Initialize preCICE internally
  precice.initialize();

  // write initial writeData to preCICE if required
  if (precice.isActionRequired(precice::constants::actionWriteInitialData()))
    {
      // store initial write_data for precice in write_data
      format_deal_to_precice(deal_to_precice);

      precice.writeBlockScalarData(write_data_id,
                                   n_interface_nodes,
                                   interface_nodes_ids.data(),
                                   write_data.data());

      precice.markActionFulfilled(precice::constants::actionWriteInitialData());

      precice.initializeData();
    }

  // read initial readData from preCICE if required for the first time step
  // FIXME: Data is already exchanged here
  if (precice.isReadDataAvailable())
    {
      precice.readBlockScalarData(read_data_id,
                                  n_interface_nodes,
                                  interface_nodes_ids.data(),
                                  read_data.data());

      format_precice_to_deal(precice_to_deal);
    }
}



template <int dim, typename VectorType, typename ParameterClass>
void
Adapter<dim, VectorType, ParameterClass>::advance(
  const VectorType &deal_to_precice,
  VectorType &      precice_to_deal,
  const double      computed_timestep_length)
{
  // This is essentially the same as during initialization
  // We have already all IDs and just need to convert our obtained data to
  // the preCICE compatible 'write_data' vector, which is done in the
  // format_deal_to_precice function. All this is of course only done in
  // case write data is required.
  if (precice.hasData(write_data_name, mesh_id))
    if (precice.isWriteDataRequired(computed_timestep_length))
      {
        format_deal_to_precice(deal_to_precice);
        precice.writeBlockScalarData(write_data_id,
                                     n_interface_nodes,
                                     interface_nodes_ids.data(),
                                     write_data.data());
      }

  // Here, we need to specify the computed time step length and pass it to
  // preCICE
  precice.advance(computed_timestep_length);

  // Here, we obtain data from another participant. Again, we insert the
  // data in our global vector by calling format_precice_to_deal

  // FIXME: isReadDataAvailable returns false, though it should be true for
  // participant b
  if (precice.hasData(read_data_name, mesh_id))
    //      if (precice.isReadDataAvailable())
    {
      precice.readBlockScalarData(read_data_id,
                                  n_interface_nodes,
                                  interface_nodes_ids.data(),
                                  read_data.data());

      format_precice_to_deal(precice_to_deal);
    }
}



template <int dim, typename VectorType, typename ParameterClass>
void
Adapter<dim, VectorType, ParameterClass>::format_deal_to_precice(
  const VectorType &deal_to_precice)
{
  // Assumption: x index is in the same position as y index in each IndexSet
  // In general, higher order support points in the element are first
  // ordered in the x component. An IndexSet for the first component might
  // look like this: [1] [3456] [11] for a 7th order 1d interface/2d cell.
  // Therefore, an index for the respective x component dof is not always
  // followed by an index on the same position for the y component

  auto dof_component = coupling_dofs.begin();
  for (int i = 0; i < n_interface_nodes; ++i)
    {
      AssertIndexRange(i, write_data.size());
      write_data[i] = deal_to_precice[*dof_component];
      ++dof_component;
    }
}



template <int dim, typename VectorType, typename ParameterClass>
void
Adapter<dim, VectorType, ParameterClass>::format_precice_to_deal(
  VectorType &precice_to_deal) const
{
  // This is the opposite direction as above. See comment there.
  auto dof_component = coupling_dofs.begin();
  for (int i = 0; i < n_interface_nodes; ++i)
    {
      AssertIndexRange(i, read_data.size());
      precice_to_deal[*dof_component] = read_data[i];
      ++dof_component;
    }
}



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
  Vector<double> dummy_vector;

  CouplingParamters                               parameters;
  const unsigned int                              interface_boundary_id;
  Adapter<dim, Vector<double>, CouplingParamters> adapter;
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
  , interface_boundary_id(3)
  , adapter(parameters, interface_boundary_id)
{}


template <int dim>
void
LaplaceProblem<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  triangulation.refine_global(3);

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
  dummy_vector.reinit(dof_handler.n_dofs());
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
  assemble_system();
  adapter.initialize(dof_handler, solution, dummy_vector);
  while (adapter.precice.isCouplingOngoing())
    {
      solve();

      adapter.advance(solution, dummy_vector, 1);
      output_results();
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  LaplaceProblem<2> laplace_problem;
  laplace_problem.run();

  return 0;
}
