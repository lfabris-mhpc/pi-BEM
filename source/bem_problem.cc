#include <deal.II/numerics/error_estimator.h>

#include <iomanip>
#include <iostream>

#include "../include/bem_problem.h"
#include "../include/constrained_matrix_complex.h"
#include "../include/laplace_kernel.h"
#include "Teuchos_TimeMonitor.hpp"

using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;

#define ENTRY EntryRaiiObject obj##LINE(__FUNCTION__);

struct EntryRaiiObject
{
  EntryRaiiObject(const char *f)
    : f_(f)
  {
    printf("Entered into %s\n", f_);
  }

  ~EntryRaiiObject()
  {
    printf("Exited from %s\n", f_);
  }

  const char *f_;
};

RCP<Time> ConstraintsTime =
  Teuchos::TimeMonitor::getNewTimer("Compute Constraints Time");
RCP<Time> AssembleTime = Teuchos::TimeMonitor::getNewTimer("Assemble Time");
RCP<Time> NormalsTime  = Teuchos::TimeMonitor::getNewTimer("Normals Time");
RCP<Time> SurfaceGradientTime =
  Teuchos::TimeMonitor::getNewTimer("SurfaceGradientTime Time");
RCP<Time> GradientTime = Teuchos::TimeMonitor::getNewTimer("Gradient Time");
RCP<Time> LacSolveTime = Teuchos::TimeMonitor::getNewTimer("LAC Solve Time");
RCP<Time> ReinitTime =
  Teuchos::TimeMonitor::getNewTimer("BEM Reinitialisation Time");

// @sect4{BEMProblem::BEMProblem and
// BEMProblem::read_parameters}
// The constructor initializes the
// variuous object in much the same
// way as done in the finite element
// programs such as step-4 or
// step-6. The only new ingredient
// here is the ParsedFunction object,
// which needs, at construction time,
// the specification of the number of
// components.
//
// For the exact solution the number
// of vector components is one, and
// no action is required since one is
// the default value for a
// ParsedFunction object. The wind,
// however, requires dim components
// to be specified. Notice that when
// declaring entries in a parameter
// file for the expression of the
// Functions::ParsedFunction, we need
// to specify the number of
// components explicitly, since the
// function
// Functions::ParsedFunction::declare_parameters
// is static, and has no knowledge of
// the number of components.
template <>
BEMProblem<3>::BEMProblem(ComputationalDomain<3> &comp_dom,
                          MPI_Comm                comm,
                          unsigned int            n_components)
  : n_components(n_components)
  , current_component(0)
  , pcout(std::cout)
  , comp_dom(comp_dom)
  , parsed_fe("Scalar FE", "FE_Q(1)")
  , parsed_gradient_fe("Vector FE", "FESystem[FE_Q(1)^3]", "u,u,u", 3)
  , dh(comp_dom.tria)
  , gradient_dh(comp_dom.tria)
  , mpi_communicator(comm)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , vector_gradients_solutions(n_components)
  , vector_surface_gradients_solutions(n_components)
{
  // Only output on first processor.
  pcout.set_condition(this_mpi_process == 0);
}

template <>
BEMProblem<2>::BEMProblem(ComputationalDomain<2> &comp_dom,
                          MPI_Comm                comm,
                          unsigned int            n_components)
  : n_components(n_components)
  , current_component(0)
  , pcout(std::cout)
  , comp_dom(comp_dom)
  , parsed_fe("Scalar FE", "FE_Q(1)")
  , parsed_gradient_fe("Vector FE", "FESystem[FE_Q(1)^2]", "u,u", 2)
  , dh(comp_dom.tria)
  , gradient_dh(comp_dom.tria)
  , mpi_communicator(comm)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , vector_gradients_solutions(n_components)
  , vector_surface_gradients_solutions(n_components)
{
  // Only output on first processor.
  pcout.set_condition(this_mpi_process == 0);
}

template <int dim>
void
BEMProblem<dim>::reinit()
{
  // ENTRY
  Teuchos::TimeMonitor LocalTimer(*ReinitTime);

  fe          = parsed_fe();
  gradient_fe = parsed_gradient_fe();
  pcout << "FE name " << fe->get_name() << std::endl;

  dh.distribute_dofs(*fe);
  gradient_dh.distribute_dofs(*gradient_fe);

  // we should choose the appropriate renumbering strategy and then stick with
  // it. in step 32 they use component_wise which is very straight-forward but
  // maybe the quickest is subdomain_wise (step 17, 18)
  DoFRenumbering::component_wise(dh);
  DoFRenumbering::component_wise(gradient_dh);

  pcout << "re-ordering vector" << std::endl;

  compute_reordering_vectors();

  DoFRenumbering::subdomain_wise(dh);
  DoFRenumbering::subdomain_wise(gradient_dh);

  vector_constraints.reinit();
  DoFTools::make_hanging_node_constraints(gradient_dh, vector_constraints);
  vector_constraints.close();
  if (mapping_type == "FE")
    {
      map_vector.reinit(gradient_dh.n_dofs());
      // Fills the euler vector with information from the Triangulation
      VectorTools::get_position_vector(gradient_dh, map_vector);
      vector_constraints.distribute(map_vector);
    }

  if (!mapping)
    {
      if (comp_dom.spheroid_bool && comp_dom.used_spherical_manifold)
        {
          for (types::global_dof_index ii = 0; ii < gradient_dh.n_dofs() / dim;
               ++ii)
            {
              map_vector[vec_original_to_sub_wise[ii]] *=
                comp_dom.spheroid_x_axis;
              map_vector[vec_original_to_sub_wise[ii + gradient_dh.n_dofs() /
                                                         dim]] *=
                comp_dom.spheroid_y_axis;
              if (dim == 3)
                {
                  map_vector[vec_original_to_sub_wise
                               [ii + gradient_dh.n_dofs() / dim]] *=
                    comp_dom.spheroid_z_axis;
                }
            }
        }

      if (mapping_type == "FE")
        {
          mapping = std::make_shared<MappingFEField<dim - 1, dim>>(gradient_dh,
                                                                   map_vector);
        }
      else
        {
          mapping = std::make_shared<MappingQ<dim - 1, dim>>(mapping_degree);
        }
    }

  pcout << "phi dofs: " << dh.n_dofs()
        << " gradient phi dofs: " << gradient_dh.n_dofs() << std::endl;
  std::vector<types::subdomain_id> dofs_domain_association(dh.n_dofs());

  DoFTools::get_subdomain_association(dh, dofs_domain_association);
  std::vector<types::subdomain_id> vector_dofs_domain_association(
    gradient_dh.n_dofs());

  DoFTools::get_subdomain_association(gradient_dh,
                                      vector_dofs_domain_association);

  this_cpu_set.clear();
  vector_this_cpu_set.clear();
  this_cpu_set.set_size(dh.n_dofs());
  vector_this_cpu_set.set_size(gradient_dh.n_dofs());

  // We compute this two vector in order to use an eventual
  // DoFRenumbering::subdomain_wise At the time being we don't. We need to
  // decide the better strategy.

  // We need to enforce consistency between the non-ghosted IndexSets.
  // To be changed accordingly with the DoFRenumbering strategy.
  pcout << "you are using " << sizeof(dh.n_dofs()) << " bytes indices"
        << std::endl;
  pcout << "setting cpu_sets" << std::endl;

  for (types::global_dof_index i = 0; i < dh.n_dofs(); ++i)
    if (dofs_domain_association[i] == this_mpi_process)
      {
        this_cpu_set.add_index(i);
        types::global_dof_index dummy = sub_wise_to_original[i];
        for (unsigned int idim = 0; idim < dim; ++idim)
          {
            vector_this_cpu_set.add_index(
              vec_original_to_sub_wise[gradient_dh.n_dofs() / dim * idim +
                                       dummy]);
          }
      }

  this_cpu_set.compress();
  vector_this_cpu_set.compress();

  // At this point we just need to create a ghosted IndexSet for the scalar
  // DoFHandler. This can be through the builtin dealii function.
  // this_cpu_set.print(std::cout);
  MPI_Barrier(mpi_communicator);
  ghosted_set.clear();
  ghosted_set.set_size(dh.n_dofs());
  ghosted_set =
    DoFTools::dof_indices_with_subdomain_association(dh, this_mpi_process);
  ghosted_set.compress();

  // standard TrilinosWrappers::MPI::Vector reinitialization.
  system_rhs.reinit(this_cpu_set, mpi_communicator);
  system_rhs_imag.reinit(this_cpu_set, mpi_communicator);
  sol.reinit(this_cpu_set, mpi_communicator);
  alpha.reinit(this_cpu_set, mpi_communicator);
  serv_phi.reinit(this_cpu_set, mpi_communicator);
  serv_phi_imag.reinit(this_cpu_set, mpi_communicator);
  serv_dphi_dn.reinit(this_cpu_set, mpi_communicator);
  serv_dphi_dn_imag.reinit(this_cpu_set, mpi_communicator);
  serv_tmp_rhs.reinit(this_cpu_set, mpi_communicator);
  serv_tmp_rhs_imag.reinit(this_cpu_set, mpi_communicator);

  // TrilinosWrappers::SparsityPattern for the BEM matricesreinitialization
  pcout << "re-initializing sparsity patterns and matrices" << std::endl;
  if (solution_method == "Direct")
    {
      full_sparsity_pattern.reinit(this_cpu_set, mpi_communicator);

      for (auto i : this_cpu_set)
        {
          for (types::global_dof_index j = 0; j < dh.n_dofs(); ++j)
            {
              full_sparsity_pattern.add(i, j);
            }
        }

      full_sparsity_pattern.compress();
      neumann_matrix.reinit(full_sparsity_pattern);
      dirichlet_matrix.reinit(full_sparsity_pattern);
    }

  pcout << "re-initialized sparsity patterns and matrices" << std::endl;
  preconditioner_band = 100;
  preconditioner_sparsity_pattern.reinit(this_cpu_set,
                                         mpi_communicator,
                                         (types::global_dof_index)
                                           preconditioner_band);
  is_preconditioner_initialized = false;

  dirichlet_nodes.reinit(this_cpu_set, mpi_communicator);
  neumann_nodes.reinit(this_cpu_set, mpi_communicator);
  robin_nodes.reinit(this_cpu_set, mpi_communicator);
  dirichlet_flags.reinit(this_cpu_set, mpi_communicator);
  neumann_flags.reinit(this_cpu_set, mpi_communicator);
  robin_flags.reinit(this_cpu_set, mpi_communicator);
  compute_dirichlet_and_neumann_dofs_vectors();
  compute_double_nodes_set();

  fma.init_fma(dh,
               double_nodes_set,
               dirichlet_nodes,
               *mapping,
               quadrature_order,
               singular_quadrature_order);

  /* TODO: unused
  // We need a TrilinosWrappers::MPI::Vector to reinit the SparsityPattern for
  // the parallel mass matrices.
  TrilinosWrappers::MPI::Vector helper(vector_this_cpu_set, mpi_communicator);
  IndexSet                      trial_index_set;
  trial_index_set.clear();
  trial_index_set =
    DoFTools::dof_indices_with_subdomain_association(gradient_dh,
                                                     this_mpi_process);
  */
  // This is the only way we could create the SparsityPattern, through the
  // Epetramap of an existing vector.
  vector_sparsity_pattern.reinit(vector_this_cpu_set,
                                 vector_this_cpu_set,
                                 mpi_communicator);
  DoFTools::make_sparsity_pattern(gradient_dh,
                                  vector_sparsity_pattern,
                                  vector_constraints,
                                  true,
                                  this_mpi_process);
  vector_sparsity_pattern.compress();
}

template <>
const Quadrature<2> &
BEMProblem<3>::get_singular_quadrature(const unsigned int index) const
{
  Assert(index < fe->dofs_per_cell, ExcIndexRange(0, fe->dofs_per_cell, index));

  static std::vector<Quadrature<2>> quadratures;
  if (quadratures.empty())
    {
      quadratures.reserve(fe->get_unit_support_points().size());
      for (const auto &unit_support_pt : fe->get_unit_support_points())
        {
          quadratures.push_back(
            QSplit<2>(QDuffy(singular_quadrature_order, 1.), unit_support_pt));
        }
    }

  return quadratures[index];
}

template <>
const Quadrature<1> &
BEMProblem<2>::get_singular_quadrature(const unsigned int index) const
{
  Assert(index < fe->dofs_per_cell, ExcIndexRange(0, fe->dofs_per_cell, index));

  static std::vector<Quadrature<1>> quadratures;
  if (quadratures.empty())
    {
      quadratures.reserve(fe->get_unit_support_points().size());
      for (const auto &unit_support_pt : fe->get_unit_support_points())
        {
          quadratures.push_back(
            QTelles<1>(singular_quadrature_order, unit_support_pt));
        }
    }

  return quadratures[index];
}

template <int dim>
void
BEMProblem<dim>::declare_parameters(ParameterHandler &prm)
{
  // In the solver section, we set
  // all SolverControl
  // parameters. The object will then
  // be fed to the GMRES solver in
  // the solve_system() function.

  prm.enter_subsection("Solver");
  SolverControl::declare_parameters(prm);
  prm.leave_subsection();

  prm.declare_entry("Preconditioner", "ILU", Patterns::Selection("ILU|AMG"));

  prm.declare_entry("Solution method",
                    "Direct",
                    Patterns::Selection("Direct|FMA"));

  prm.enter_subsection("Quadrature rules");
  {
    prm.declare_entry("Quadrature type",
                      "gauss",
                      Patterns::Selection(
                        QuadratureSelector<(dim - 1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "4", Patterns::Integer());
    prm.declare_entry("Singular quadrature order", "5", Patterns::Integer());
  }
  prm.leave_subsection();

  prm.declare_entry("Mapping Type", "FE", Patterns::Selection("FE|Q"));

  prm.declare_entry("Mapping Q Degree", "1", Patterns::Integer());

  prm.declare_entry("Continuos gradient across edges",
                    "true",
                    Patterns::Bool());
}

template <int dim>
void
BEMProblem<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Solver");
  solver_control.parse_parameters(prm);
  prm.leave_subsection();

  preconditioner_type = prm.get("Preconditioner");

  solution_method = prm.get("Solution method");

  prm.enter_subsection("Quadrature rules");
  {
    quadrature = std::shared_ptr<Quadrature<dim - 1>>(
      new QuadratureSelector<dim - 1>(prm.get("Quadrature type"),
                                      prm.get_integer("Quadrature order")));
    quadrature_order          = prm.get_integer("Quadrature order");
    singular_quadrature_order = prm.get_integer("Singular quadrature order");
  }
  prm.leave_subsection();

  mapping_type       = prm.get("Mapping Type");
  mapping_degree     = prm.get_integer("Mapping Q Degree");
  continuos_gradient = prm.get_bool("Continuos gradient across edges");
}

template <int dim>
void
BEMProblem<dim>::compute_dirichlet_and_neumann_dofs_vectors()
{
  can_determine_phi = false;

  Vector<double> non_partitioned_dirichlet_nodes(dh.n_dofs());
  Vector<double> non_partitioned_neumann_nodes(dh.n_dofs());
  Vector<double> non_partitioned_robin_nodes(dh.n_dofs());

  Vector<double> non_partitioned_dirichlet_flags(dh.n_dofs());
  Vector<double> non_partitioned_neumann_flags(dh.n_dofs());
  Vector<double> non_partitioned_robin_flags(dh.n_dofs());

  // defaulting to neumann
  // vector_shift(non_partitioned_neumann_nodes, 1.);
  std::vector<types::global_dof_index> dofs(fe->dofs_per_cell);
  unsigned int                         local_can_determine_phi = 0;

  for (const auto &cell : dh.active_cell_iterators())
    {
      if (cell->subdomain_id() == this_mpi_process)
        {
          bool is_dirichlet = std::find(comp_dom.dirichlet_boundary_ids.begin(),
                                        comp_dom.dirichlet_boundary_ids.end(),
                                        cell->boundary_id()) !=
                              comp_dom.dirichlet_boundary_ids.end();
          if (is_dirichlet)
            {
              cell->get_dof_indices(dofs);
              for (auto i : dofs)
                {
                  non_partitioned_dirichlet_flags(i) = 1;

                  // mark dofs on masking vectors
                  non_partitioned_dirichlet_nodes(i) = 1;
                  non_partitioned_neumann_nodes(i)   = 0;
                  non_partitioned_robin_nodes(i)     = 0;
                }

              local_can_determine_phi = 1;
            }
          else
            {
              bool is_neumann = std::find(comp_dom.neumann_boundary_ids.begin(),
                                          comp_dom.neumann_boundary_ids.end(),
                                          cell->boundary_id()) !=
                                comp_dom.neumann_boundary_ids.end();
              if (is_neumann)
                {
                  cell->get_dof_indices(dofs);
                  for (auto i : dofs)
                    {
                      non_partitioned_neumann_flags(i) = 1;

                      // mark dofs on masking vectors
                      non_partitioned_dirichlet_nodes(i) = 0;
                      non_partitioned_neumann_nodes(i)   = 1;
                      non_partitioned_robin_nodes(i)     = 0;
                    }
                }
              else
                {
#ifdef DEBUG
                  bool is_robin = std::find(comp_dom.robin_boundary_ids.begin(),
                                            comp_dom.robin_boundary_ids.end(),
                                            cell->boundary_id()) !=
                                  comp_dom.robin_boundary_ids.end();
                  Assert(is_robin, ExcInternalError());
#endif
                  cell->get_dof_indices(dofs);
                  for (auto i : dofs)
                    {
                      non_partitioned_robin_flags(i) = 1;

                      // mark dofs on masking vectors
                      if (!non_partitioned_dirichlet_nodes(i) &&
                          !non_partitioned_neumann_nodes(i))
                        {
                          non_partitioned_dirichlet_nodes(i) = 0;
                          non_partitioned_neumann_nodes(i)   = 0;
                          non_partitioned_robin_nodes(i)     = 1;
                        }

                      local_can_determine_phi = 1;
                    }
                }
            }
        }
    }

  for (auto i : this_cpu_set)
    {
      dirichlet_nodes(i) = non_partitioned_dirichlet_nodes(i);
      neumann_nodes(i)   = non_partitioned_neumann_nodes(i);
      robin_nodes(i)     = non_partitioned_robin_nodes(i);

      dirichlet_flags(i) = non_partitioned_dirichlet_flags(i);
      neumann_flags(i)   = non_partitioned_neumann_flags(i);
      robin_flags(i)     = non_partitioned_robin_flags(i);
    }

  {
    Vector<double> localized_dirichlet(dirichlet_nodes);
    pcout << "Number of Dirichlet dofs: "
          << (int)(localized_dirichlet.size() *
                   localized_dirichlet.mean_value())
          << std::endl;
    Vector<double> localized_neumann(neumann_nodes);
    pcout << "Number of Neumann dofs: "
          << (int)(localized_neumann.size() * localized_neumann.mean_value())
          << std::endl;
    Vector<double> localized_robin(robin_nodes);
    pcout << "Number of Robin dofs: "
          << (int)(localized_robin.size() * localized_robin.mean_value())
          << std::endl;

    Vector<double> localized_dirichlet2(dirichlet_flags);
    pcout << "Number of Dirichlet flags: "
          << (int)(localized_dirichlet2.size() *
                   localized_dirichlet2.mean_value())
          << std::endl;
    Vector<double> localized_neumann2(neumann_flags);
    pcout << "Number of Neumann flags: "
          << (int)(localized_neumann2.size() * localized_neumann2.mean_value())
          << std::endl;
    Vector<double> localized_robin2(robin_flags);
    pcout << "Number of Robin flags: "
          << (int)(localized_robin2.size() * localized_robin2.mean_value())
          << std::endl;
  }

  unsigned int global_can_determine_phi;
  MPI_Allreduce(&local_can_determine_phi,
                &global_can_determine_phi,
                1,
                MPI_UNSIGNED,
                MPI_MAX,
                mpi_communicator);

  if (global_can_determine_phi > 0)
    {
      can_determine_phi = true;
    }
}

template <int dim>
void
BEMProblem<dim>::compute_double_nodes_set()
{
  double tol = 1e-10;
  double_nodes_set.clear();
  double_nodes_set.resize(dh.n_dofs());
  std::vector<Point<dim>> support_points(dh.n_dofs());

  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*mapping,
                                                     dh,
                                                     support_points);
  std::vector<types::global_dof_index> face_dofs(fe->dofs_per_face);

  edge_set.clear();
  edge_set.set_size(dh.n_dofs());
  for (const auto &cell : dh.active_cell_iterators())
    {
      for (unsigned int f = 0; f < GeometryInfo<dim - 1>::faces_per_cell; ++f)
        {
          if (cell->face(f)->at_boundary())
            {
              cell->face(f)->get_dof_indices(face_dofs);
              edge_set.add_indices(face_dofs.cbegin(), face_dofs.cend());
            }
        }
    }
  edge_set.compress();

  for (types::global_dof_index i = 0; i < dh.n_dofs(); ++i)
    {
      double_nodes_set[i].insert(i);
    }
  for (auto i : edge_set)
    {
      auto jiter = edge_set.at(i);
      ++jiter;
      for (; jiter != edge_set.end(); ++jiter)
        {
          const auto j = *jiter;
          if (support_points[i].distance_square(support_points[j]) <
              (tol * tol))
            {
              double_nodes_set[i].insert(j);
              double_nodes_set[j].insert(i);
            }
        }
      /*
      for (auto j : edge_set)
        {
          if (support_points[i].distance_square(support_points[j]) <
              (tol * tol))
            {
              double_nodes_set[i].insert(j);
            }
        }
      */
    }
}

template <int dim>
void
BEMProblem<dim>::compute_reordering_vectors()
{
  original_to_sub_wise.resize(dh.n_dofs());
  sub_wise_to_original.resize(dh.n_dofs());
  vec_original_to_sub_wise.resize(gradient_dh.n_dofs());
  vec_sub_wise_to_original.resize(gradient_dh.n_dofs());

  DoFRenumbering::compute_subdomain_wise(original_to_sub_wise, dh);
  DoFRenumbering::compute_subdomain_wise(vec_original_to_sub_wise, gradient_dh);

  for (types::global_dof_index i = 0; i < gradient_dh.n_dofs(); ++i)
    {
      if (i < dh.n_dofs())
        {
          sub_wise_to_original[original_to_sub_wise[i]] = i;
        }
      vec_sub_wise_to_original[vec_original_to_sub_wise[i]] = i;
    }
}

template <int dim>
void
BEMProblem<dim>::assemble_system()
{
  Teuchos::TimeMonitor LocalTimer(*AssembleTime);
  pcout << "(Directly) Assembling system matrices" << std::endl;

  neumann_matrix   = 0;
  dirichlet_matrix = 0;

  // Next, we initialize an FEValues
  // object with the quadrature
  // formula for the integration of
  // the kernel in non singular
  // cells. This quadrature is
  // selected with the parameter
  // file, and needs to be quite
  // precise, since the functions we
  // are integrating are not
  // polynomial functions.
  FEValues<dim - 1, dim> fe_v(*mapping,
                              *fe,
                              *quadrature,
                              update_values | update_normal_vectors |
                                update_quadrature_points | update_JxW_values);

  const unsigned int                   n_q_points = fe_v.n_quadrature_points;
  std::vector<types::global_dof_index> local_dof_indices(fe->dofs_per_cell);
  pcout << "DoFs per cell: " << fe->dofs_per_cell << " " << std::endl;

  // Unlike in finite element
  // methods, if we use a collocation
  // boundary element method, then in
  // each assembly loop we only
  // assemble the information that
  // refers to the coupling between
  // one degree of freedom (the
  // degree associated with support
  // point $i$) and the current
  // cell. This is done using a
  // vector of fe->dofs_per_cell
  // elements, which will then be
  // distributed to the matrix in the
  // global row $i$. The following
  // object will hold this
  // information:
  Vector<double> local_neumann_matrix_row_i(fe->dofs_per_cell);
  Vector<double> local_dirichlet_matrix_row_i(fe->dofs_per_cell);

  // Now that we have checked that
  // the number of vertices is equal
  // to the number of degrees of
  // freedom, we construct a vector
  // of support points which will be
  // used in the local integrations:
  std::vector<Point<dim>> support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*mapping,
                                                     dh,
                                                     support_points);

  // After doing so, we can start the
  // integration loop over all cells,
  // where we first initialize the
  // FEValues object and get the
  // values of $\mathbf{\tilde v}$ at
  // the quadrature points (this
  // vector field should be constant,
  // but it doesn't hurt to be more
  // general):

  Point<dim> D;
  double     s;

  for (const auto &cell : dh.active_cell_iterators())
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      // std::vector<Point> and std::vector<Tensor>
      const auto &q_points = fe_v.get_quadrature_points();
      const auto &normals  = fe_v.get_normal_vectors();

      // We then form the integral over
      // the current cell for all
      // degrees of freedom (note that
      // this includes degrees of
      // freedom not located on the
      // current cell, a deviation from
      // the usual finite element
      // integrals). The integral that
      // we need to perform is singular
      // if one of the local degrees of
      // freedom is the same as the
      // support point $i$. A the
      // beginning of the loop we
      // therefore check wether this is
      // the case, and we store which
      // one is the singular index:
      for (auto i : this_cpu_set)
        { // these must now be the locally owned dofs.
          // the rest should stay the same
          local_neumann_matrix_row_i   = 0;
          local_dirichlet_matrix_row_i = 0;

          bool         is_singular    = false;
          unsigned int singular_index = numbers::invalid_unsigned_int;

          // is any dof of the current cell, a duplicate of i?
          for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
            {
              if (double_nodes_set[i].count(local_dof_indices[j]) > 0)
                {
                  singular_index = j;
                  is_singular    = true;
                  break;
                }
            }

          // We then perform the
          // integral. If the index $i$
          // is not one of the local
          // degrees of freedom, we
          // simply have to add the
          // single layer terms to the
          // right hand side, and the
          // double layer terms to the
          // matrix:
          if (!is_singular)
            {
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  // const Tensor<1, dim> R = q_points[q] - support_points[i];
                  LaplaceKernel::kernels(q_points[q] - support_points[i], D, s);

                  for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
                    {
                      const auto tmp = fe_v.shape_value(j, q) * fe_v.JxW(q);

                      local_neumann_matrix_row_i(j) += ((D * normals[q]) * tmp);
                      local_dirichlet_matrix_row_i(j) += (s * tmp);
                    }
                }
            }
          else
            {
              // Now we treat the more
              // delicate case. If we
              // are here, this means
              // that the cell that
              // runs on the $j$ index
              // contains
              // support_point[i]. In
              // this case both the
              // single and the double
              // layer potential are
              // singular, and they
              // require special
              // treatment.
              //
              // Whenever the
              // integration is
              // performed with the
              // singularity inside the
              // given cell, then a
              // special quadrature
              // formula is used that
              // allows one to
              // integrate arbitrary
              // functions against a
              // singular weight on the
              // reference cell.
              // Notice that singular
              // integration requires a
              // careful selection of
              // the quadrature
              // rules. In particular
              // the deal.II library
              // provides quadrature
              // rules which are
              // taylored for
              // logarithmic
              // singularities
              // (QGaussLog,
              // QGaussLogR), as well
              // as for 1/R
              // singularities
              // (QGaussOneOverR).
              //
              // Singular integration
              // is typically obtained
              // by constructing
              // weighted quadrature
              // formulas with singular
              // weights, so that it is
              // possible to write
              //
              // \f[
              //   \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i)
              // \f]
              //
              // where $s(x)$ is a given
              // singularity, and the weights
              // and quadrature points
              // $w_i,q_i$ are carefully
              // selected to make the formula
              // above an equality for a
              // certain class of functions
              // $f(x)$.
              //
              // In all the finite
              // element examples we
              // have seen so far, the
              // weight of the
              // quadrature itself
              // (namely, the function
              // $s(x)$), was always
              // constantly equal to 1.
              // For singular
              // integration, we have
              // two choices: we can
              // use the definition
              // above, factoring out
              // the singularity from
              // the integrand (i.e.,
              // integrating $f(x)$
              // with the special
              // quadrature rule), or
              // we can ask the
              // quadrature rule to
              // "normalize" the
              // weights $w_i$ with
              // $s(q_i)$:
              //
              // \f[
              //   \int_K f(x) s(x) dx =
              //   \int_K g(x) dx = \sum_{i=1}^N \frac{w_i}{s(q_i)} g(q_i)
              // \f]
              //
              // We use this second
              // option, through the @p
              // factor_out_singularity
              // parameter of both
              // QGaussLogR and
              // QGaussOneOverR.
              //
              // These integrals are
              // somewhat delicate,
              // especially in two
              // dimensions, due to the
              // transformation from
              // the real to the
              // reference cell, where
              // the variable of
              // integration is scaled
              // with the determinant
              // of the transformation.
              //
              // In two dimensions this
              // process does not
              // result only in a
              // factor appearing as a
              // constant factor on the
              // entire integral, but
              // also on an additional
              // integral alltogether
              // that needs to be
              // evaluated:
              //
              // \f[
              //  \int_0^1 f(x)\ln(x/\alpha) dx =
              //  \int_0^1 f(x)\ln(x) dx - \int_0^1 f(x) \ln(\alpha) dx.
              // \f]
              //
              // This process is taken care of by
              // the constructor of the QGaussLogR
              // class, which adds additional
              // quadrature points and weights to
              // take into consideration also the
              // second part of the integral.
              //
              // A similar reasoning
              // should be done in the
              // three dimensional
              // case, since the
              // singular quadrature is
              // taylored on the
              // inverse of the radius
              // $r$ in the reference
              // cell, while our
              // singular function
              // lives in real space,
              // however in the three
              // dimensional case
              // everything is simpler
              // because the
              // singularity scales
              // linearly with the
              // determinant of the
              // transformation. This
              // allows us to build the
              // singular two
              // dimensional quadrature
              // rules once and for all
              // outside the loop over
              // all cells, using only
              // a pointer where needed.
              //
              // Notice that in one
              // dimensional
              // integration this is
              // not possible, since we
              // need to know the
              // scaling parameter for
              // the quadrature, which
              // is not known a
              // priori. Here, the
              // quadrature rule itself
              // depends also on the
              // size of the current
              // cell. For this reason,
              // it is necessary to
              // create a new
              // quadrature for each
              // singular
              // integration. Since we
              // create it using the
              // new operator of C++,
              // we also need to
              // destroy it using the
              // dual of new:
              // delete. This is done
              // at the end, and only
              // if dim == 2.
              //
              // Putting all this into a
              // dimension independent
              // framework requires a little
              // trick. The problem is that,
              // depending on dimension, we'd
              // like to either assign a
              // QGaussLogR<1> or a
              // QGaussOneOverR<2> to a
              // Quadrature<dim-1>. C++
              // doesn't allow this right
              // away, and neither is a
              // static_cast
              // possible. However, we can
              // attempt a dynamic_cast: the
              // implementation will then
              // look up at run time whether
              // the conversion is possible
              // (which we <em>know</em> it
              // is) and if that isn't the
              // case simply return a null
              // pointer. To be sure we can
              // then add a safety check at
              // the end:
              Assert(singular_index != numbers::invalid_unsigned_int,
                     ExcInternalError());

              // pointer trick
              const Quadrature<dim - 1> *singular_quadrature =
                &(get_singular_quadrature(singular_index));
              Assert(singular_quadrature, ExcInternalError());

              FEValues<dim - 1, dim> fe_v_singular(*mapping,
                                                   *fe,
                                                   *singular_quadrature,
                                                   update_jacobians |
                                                     update_values |
                                                     update_normal_vectors |
                                                     update_quadrature_points);
              fe_v_singular.reinit(cell);

              // std::vector<Point> and std::vector<Tensor>
              const auto &singular_normals = fe_v_singular.get_normal_vectors();
              const auto &singular_q_points =
                fe_v_singular.get_quadrature_points();

              for (unsigned int q = 0; q < singular_quadrature->size(); ++q)
                {
                  // const Tensor<1, dim> R = singular_q_points[q] -
                  // support_points[i];
                  LaplaceKernel::kernels(
                    singular_q_points[q] - support_points[i], D, s);

                  for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
                    {
                      const auto tmp =
                        fe_v_singular.shape_value(j, q) * fe_v_singular.JxW(q);

                      local_neumann_matrix_row_i(j) +=
                        ((D * singular_normals[q]) * tmp);
                      local_dirichlet_matrix_row_i(j) += (s * tmp);
                    }
                }
            }

          // Finally, we need to add
          // the contributions of the
          // current cell to the
          // global matrix.
          for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
            {
              neumann_matrix.add(i,
                                 local_dof_indices[j],
                                 local_neumann_matrix_row_i(j));
              dirichlet_matrix.add(i,
                                   local_dof_indices[j],
                                   local_dirichlet_matrix_row_i(j));
            }
        }
    }

  // The second part of the integral
  // operator is the term
  // $\alpha(\mathbf{x}_i)
  // \phi_j(\mathbf{x}_i)$. Since we
  // use a collocation scheme,
  // $\phi_j(\mathbf{x}_i)=\delta_{ij}$
  // and the corresponding matrix is
  // a diagonal one with entries
  // equal to $\alpha(\mathbf{x}_i)$.

  // One quick way to compute this
  // diagonal matrix of the solid
  // angles, is to use the Neumann
  // matrix itself. It is enough to
  // multiply the matrix with a
  // vector of elements all equal to
  // -1, to get the diagonal matrix
  // of the alpha angles, or solid
  // angles (see the formula in the
  // introduction for this). The
  // result is then added back onto
  // the system matrix object to
  // yield the final form of the
  // matrix:

  pcout << "done assembling system matrices" << std::endl;
}

template <int dim>
void
BEMProblem<dim>::compute_alpha()
{
  static TrilinosWrappers::MPI::Vector ones, zeros, dummy;
  if (ones.size() != dh.n_dofs())
    {
      ones.reinit(this_cpu_set, mpi_communicator);
      vector_shift(ones, -1.);
      zeros.reinit(this_cpu_set, mpi_communicator);
      dummy.reinit(this_cpu_set, mpi_communicator);
    }

  if (solution_method == "Direct")
    {
      neumann_matrix.vmult(alpha, ones);
    }
  else
    {
      AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));

      fma.generate_multipole_expansions(ones, zeros);
      fma.multipole_matr_vect_products(ones, zeros, alpha, dummy);
    }
}

template <int dim>
void
BEMProblem<dim>::vmult(TrilinosWrappers::MPI::Vector &      dst,
                       const TrilinosWrappers::MPI::Vector &src) const
{
  // the Robin nodes participate with their unknowns carrying phi
  //(alpha + N) * serv_phi - D * serv_dphi_dn
  // becomes
  //(alpha + N) * (serv_phi + serv_phi_robin) - D * (serv_dphi_dn -
  // robin_matrix_diagonal.scale(serv_phi_robin))

  serv_phi = src;
  if (!can_determine_phi)
    {
      vector_shift(serv_phi, -serv_phi.l2_norm());
    }
  serv_dphi_dn   = src;
  serv_phi_robin = serv_phi;

  TrilinosWrappers::MPI::Vector matrVectProdN;
  TrilinosWrappers::MPI::Vector matrVectProdD;

  matrVectProdN.reinit(this_cpu_set, mpi_communicator);
  matrVectProdD.reinit(this_cpu_set, mpi_communicator);

  dst = 0;

  serv_phi.scale(neumann_nodes);
  serv_dphi_dn.scale(dirichlet_nodes);
  serv_phi_robin.scale(robin_nodes);

  if (solution_method == "Direct")
    {
      serv_phi += serv_phi_robin;
      serv_phi_robin.scale(robin_matrix_diagonal);
      serv_dphi_dn -= serv_phi_robin;

      dirichlet_matrix.vmult(dst, serv_dphi_dn);
      dst *= -1;
      neumann_matrix.vmult_add(dst, serv_phi);
      serv_phi.scale(alpha);
      dst += serv_phi;
    }
  else
    {
      AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));

      fma.generate_multipole_expansions(serv_phi, serv_dphi_dn);

      serv_phi += serv_phi_robin;
      serv_phi_robin.scale(robin_matrix_diagonal);
      serv_dphi_dn -= serv_phi_robin;

      fma.multipole_matr_vect_products(serv_phi,
                                       serv_dphi_dn,
                                       matrVectProdN,
                                       matrVectProdD);
      dst += matrVectProdD;
      dst *= -1;
      dst += matrVectProdN;
      serv_phi.scale(alpha);
      dst += serv_phi;
    }

  if (!can_determine_phi)
    {
      vector_shift(dst, -dst.l2_norm());
    }
  dst.compress(VectorOperation::add);
}

template <int dim>
void
BEMProblem<dim>::vmult(TrilinosWrappers::MPI::Vector &      dst,
                       TrilinosWrappers::MPI::Vector &      dst_imag,
                       const TrilinosWrappers::MPI::Vector &src,
                       const TrilinosWrappers::MPI::Vector &src_imag) const
{
  // the Robin nodes participate with their unknowns carrying phi
  //(alpha + N) * serv_phi - D * serv_dphi_dn
  // becomes
  //(alpha + N) * (serv_phi + serv_phi_robin) - D * (serv_dphi_dn -
  // robin_matrix_diagonal.scale(serv_phi_robin))

  serv_phi      = src;
  serv_phi_imag = src_imag;
  if (!can_determine_phi)
    {
      auto shift = std::sqrt(serv_phi.norm_sqr() + serv_phi_imag.l2_norm());
      vector_shift(serv_phi, -shift);
      vector_shift(serv_phi_imag, -shift);
    }
  serv_dphi_dn        = src;
  serv_dphi_dn_imag   = src_imag;
  serv_phi_robin      = serv_phi;
  serv_phi_robin_imag = serv_phi_imag;

  TrilinosWrappers::MPI::Vector matrVectProdN;
  TrilinosWrappers::MPI::Vector matrVectProdN_imag;
  TrilinosWrappers::MPI::Vector matrVectProdD;
  TrilinosWrappers::MPI::Vector matrVectProdD_imag;
  TrilinosWrappers::MPI::Vector tmp1, tmp2;

  matrVectProdN.reinit(this_cpu_set, mpi_communicator);
  matrVectProdN_imag.reinit(this_cpu_set, mpi_communicator);
  matrVectProdD.reinit(this_cpu_set, mpi_communicator);
  matrVectProdD_imag.reinit(this_cpu_set, mpi_communicator);
  tmp1.reinit(this_cpu_set, mpi_communicator);
  tmp2.reinit(this_cpu_set, mpi_communicator);

  dst      = 0;
  dst_imag = 0;

  serv_phi.scale(neumann_nodes);
  serv_phi_imag.scale(neumann_nodes);
  serv_dphi_dn_imag.scale(dirichlet_nodes);
  serv_dphi_dn.scale(dirichlet_nodes);
  serv_phi_robin.scale(robin_nodes);
  serv_phi_robin_imag.scale(robin_nodes);

  if (solution_method == "Direct")
    {
      serv_phi += serv_phi_robin;
      serv_phi_imag += serv_phi_robin_imag;
      // robin_matrix_diagonal is complex
      // serv_dphi_dn += -robin_matrix_diagonal*serv_phi +
      // robin_matrix_diagonal_imag * serv_phi_imag
      // and
      // serv_dphi_dn_imag += -robin_matrix_diagonal_imag*serv_phi -
      // robin_matrix_diagonal * serv_phi_imag
      tmp1 = serv_phi_robin;
      tmp1.scale(robin_matrix_diagonal);
      tmp2 = serv_phi_robin_imag;
      tmp2.scale(robin_matrix_diagonal_imag);
      serv_dphi_dn -= tmp1;
      serv_dphi_dn += tmp2;

      // these can be destructive
      serv_phi_robin.scale(robin_matrix_diagonal_imag);
      serv_phi_robin_imag.scale(robin_matrix_diagonal);
      serv_dphi_dn_imag -= serv_phi_robin;
      serv_dphi_dn_imag -= serv_phi_robin_imag;

      dirichlet_matrix.vmult(dst, serv_dphi_dn);
      dirichlet_matrix.vmult(dst_imag, serv_dphi_dn_imag);
      dst *= -1;
      dst_imag *= -1;
      neumann_matrix.vmult_add(dst, serv_phi);
      neumann_matrix.vmult_add(dst_imag, serv_phi_imag);
      serv_phi.scale(alpha);
      serv_phi_imag.scale(alpha);
      dst += serv_phi;
      dst_imag += serv_phi_imag;
    }
  else
    {
      AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));
      // TODO: this is much more involved, expansions are to be redone for each
      // product?
      fma.generate_multipole_expansions(serv_phi, serv_dphi_dn);
      serv_phi += serv_phi_robin;
      serv_phi_robin.scale(robin_matrix_diagonal);
      serv_dphi_dn -= serv_phi_robin;
      fma.multipole_matr_vect_products(serv_phi,
                                       serv_dphi_dn,
                                       matrVectProdN,
                                       matrVectProdD);
      dst += matrVectProdD;
      dst *= -1;
      dst += matrVectProdN;
      serv_phi.scale(alpha);
      dst += serv_phi;
    }

  if (!can_determine_phi)
    {
      auto shift = std::sqrt(dst.norm_sqr() + dst_imag.l2_norm());
      vector_shift(dst, -shift);
      vector_shift(dst_imag, -shift);
    }
  dst.compress(VectorOperation::add);
  dst_imag.compress(VectorOperation::add);
}

template <int dim>
void
BEMProblem<dim>::compute_rhs(TrilinosWrappers::MPI::Vector &      dst,
                             const TrilinosWrappers::MPI::Vector &src) const
{
  // the Robin nodes participate with their unknowns carrying phi; the
  // inhomogeneity is accounted for as if included in dphi_dn
  //-(alpha + N) * serv_phi + D * serv_dphi_dn
  // becomes
  //-(alpha + N) * serv_phi + D * (serv_dphi_dn + robin_rhs)
  serv_phi     = src;
  serv_dphi_dn = src;

  static TrilinosWrappers::MPI::Vector matrVectProdN;
  static TrilinosWrappers::MPI::Vector matrVectProdD;

  matrVectProdN.reinit(this_cpu_set, mpi_communicator);
  matrVectProdD.reinit(this_cpu_set, mpi_communicator);

  // Robin nodes are accounted for by robin_rhs
  serv_phi.scale(dirichlet_nodes);
  serv_dphi_dn.scale(neumann_nodes);
  // cut the robin_rhs to only the true robin nodes
  serv_phi_robin = robin_rhs;
  serv_phi_robin.scale(robin_nodes);

  if (solution_method == "Direct")
    {
      neumann_matrix.vmult(dst, serv_phi);
      serv_phi.scale(alpha);
      dst += serv_phi;
      dst *= -1;
      serv_dphi_dn += serv_phi_robin;
      dirichlet_matrix.vmult_add(dst, serv_dphi_dn);
    }
  else
    {
      AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));

      fma.generate_multipole_expansions(serv_phi, serv_dphi_dn);
      serv_dphi_dn += serv_phi_robin;
      fma.multipole_matr_vect_products(serv_phi,
                                       serv_dphi_dn,
                                       matrVectProdN,
                                       matrVectProdD);
      serv_phi.scale(alpha);
      dst += matrVectProdN;
      dst += serv_phi;
      dst *= -1;
      dst += matrVectProdD;
    }
}

template <int dim>
void
BEMProblem<dim>::compute_rhs(
  TrilinosWrappers::MPI::Vector &      dst,
  TrilinosWrappers::MPI::Vector &      dst_imag,
  const TrilinosWrappers::MPI::Vector &src,
  const TrilinosWrappers::MPI::Vector &src_imag) const
{
  // the Robin nodes participate with their unknowns carrying phi; the
  // inhomogeneity is accounted for as if included in dphi_dn
  //-(alpha + N) * serv_phi + D * serv_dphi_dn
  // becomes
  //-(alpha + N) * serv_phi + D * (serv_dphi_dn + robin_rhs)
  serv_phi          = src;
  serv_phi_imag     = src_imag;
  serv_dphi_dn      = src;
  serv_dphi_dn_imag = src_imag;

  static TrilinosWrappers::MPI::Vector matrVectProdN;
  static TrilinosWrappers::MPI::Vector matrVectProdN_imag;
  static TrilinosWrappers::MPI::Vector matrVectProdD;
  static TrilinosWrappers::MPI::Vector matrVectProdD_imag;

  matrVectProdN.reinit(this_cpu_set, mpi_communicator);
  matrVectProdN_imag.reinit(this_cpu_set, mpi_communicator);
  matrVectProdD.reinit(this_cpu_set, mpi_communicator);
  matrVectProdD_imag.reinit(this_cpu_set, mpi_communicator);

  // Robin nodes are accounted for by robin_rhs
  serv_phi.scale(dirichlet_nodes);
  serv_phi_imag.scale(dirichlet_nodes);
  serv_dphi_dn.scale(neumann_nodes);
  serv_dphi_dn_imag.scale(neumann_nodes);
  // cut the robin_rhs to only the true robin nodes
  serv_phi_robin      = robin_rhs;
  serv_phi_robin_imag = robin_rhs_imag;
  serv_phi_robin.scale(robin_nodes);
  serv_phi_robin_imag.scale(robin_nodes);

  if (solution_method == "Direct")
    {
      neumann_matrix.vmult(dst, serv_phi);
      neumann_matrix.vmult(dst_imag, serv_phi_imag);
      serv_phi.scale(alpha);
      serv_phi_imag.scale(alpha);
      dst += serv_phi;
      dst_imag += serv_phi_imag;
      dst *= -1;
      dst_imag *= -1;
      serv_dphi_dn += serv_phi_robin;
      serv_dphi_dn_imag += serv_phi_robin_imag;
      dirichlet_matrix.vmult_add(dst, serv_dphi_dn);
      dirichlet_matrix.vmult_add(dst_imag, serv_dphi_dn_imag);
    }
  else
    {
      AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));

      fma.generate_multipole_expansions(serv_phi, serv_dphi_dn);
      serv_dphi_dn += serv_phi_robin;
      serv_dphi_dn_imag += serv_phi_robin_imag;
      fma.multipole_matr_vect_products(serv_phi,
                                       serv_dphi_dn,
                                       matrVectProdN,
                                       matrVectProdD);
      fma.multipole_matr_vect_products(serv_phi_imag,
                                       serv_dphi_dn_imag,
                                       matrVectProdN_imag,
                                       matrVectProdD_imag);
      serv_phi.scale(alpha);
      serv_phi_imag.scale(alpha);
      dst += matrVectProdN;
      dst_imag += matrVectProdN_imag;
      dst += serv_phi;
      dst_imag += serv_phi_imag;
      dst *= -1;
      dst_imag *= -1;
      dst += matrVectProdD;
      dst_imag += matrVectProdD_imag;
    }
}

// @sect4{BEMProblem::solve_system}

// The next function simply solves
// the linear system.
template <int dim>
void
BEMProblem<dim>::solve_system(TrilinosWrappers::MPI::Vector &      phi,
                              TrilinosWrappers::MPI::Vector &      dphi_dn,
                              const TrilinosWrappers::MPI::Vector &tmp_rhs)
{
  Teuchos::TimeMonitor                       LocalTimer(*LacSolveTime);
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(
    solver_control,
    SolverGMRES<TrilinosWrappers::MPI::Vector>::AdditionalData(100));

  system_rhs = 0;
  sol        = 0;
  alpha      = 0;

  compute_alpha();
  compute_rhs(system_rhs, tmp_rhs);

  compute_constraints(constr_cpu_set, constraints, tmp_rhs);
  ConstrainedOperator<TrilinosWrappers::MPI::Vector, BEMProblem<dim>> cc(
    *this, constraints, constr_cpu_set, mpi_communicator);

  cc.distribute_rhs(system_rhs);
  system_rhs.compress(VectorOperation::insert);

  if (solution_method == "Direct")
    {
      assemble_preconditioner();
      sol.sadd(1., 0., system_rhs);
      solver.solve(cc, sol, system_rhs, preconditioner);
    }
  else
    {
      AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));

      TrilinosWrappers::PreconditionILU &fma_preconditioner =
        fma.FMA_preconditioner(alpha, constraints);
      solver.solve(cc, sol, system_rhs, fma_preconditioner);
    }

  for (auto i : this_cpu_set)
    {
      if (neumann_nodes(i) == 1)
        {
          phi(i) = sol(i);
        }
      else if (dirichlet_nodes(i) == 1)
        {
          dphi_dn(i) = sol(i);
        }
      else
        {
          Assert(robin_nodes(i) == 1, ExcInternalError());
          phi(i)     = sol(i);
          dphi_dn(i) = robin_rhs(i) - robin_matrix_diagonal(i) * phi(i);
        }
    }

  phi(this_cpu_set.nth_index_in_set(0)) = phi(this_cpu_set.nth_index_in_set(0));
  dphi_dn(this_cpu_set.nth_index_in_set(0)) =
    dphi_dn(this_cpu_set.nth_index_in_set(0));
  phi.compress(VectorOperation::insert);
  dphi_dn.compress(VectorOperation::insert);
}

template <int dim>
void
BEMProblem<dim>::solve_system(TrilinosWrappers::MPI::Vector &      phi,
                              TrilinosWrappers::MPI::Vector &      phi_imag,
                              TrilinosWrappers::MPI::Vector &      dphi_dn,
                              TrilinosWrappers::MPI::Vector &      dphi_dn_imag,
                              const TrilinosWrappers::MPI::Vector &tmp_rhs,
                              const TrilinosWrappers::MPI::Vector &tmp_rhs_imag)
{
  Teuchos::TimeMonitor                       LocalTimer(*LacSolveTime);
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(
    solver_control,
    SolverGMRES<TrilinosWrappers::MPI::Vector>::AdditionalData(100));

  // TODO: needs to assemble and pass double-length sol and rhs
  // the ConstrainedComplexMatrix will need to split and reconstruct the
  // subvectors at each iteration
  // must take care of the subvector constraints
  // imaginary part has the same indexes as the real one, shifted by the
  // original size

  system_rhs      = 0;
  system_rhs_imag = 0;
  sol             = 0;
  alpha           = 0;

  compute_alpha();
  compute_rhs(system_rhs, system_rhs_imag, tmp_rhs, tmp_rhs_imag);

  compute_constraints(constr_cpu_set, constraints, tmp_rhs);
  set_current_phi_component(current_component + 1);
  compute_constraints(constr_cpu_set, constraints_imag, tmp_rhs_imag);
  set_current_phi_component(current_component - 1);

  // assemble the complex datastructs
  IndexSet complex_cpu_set(2 * this_cpu_set.size());
  complex_cpu_set.add_indices(this_cpu_set);
  complex_cpu_set.add_indices(this_cpu_set, this_cpu_set.size());
  complex_cpu_set.compress();
  TrilinosWrappers::MPI::Vector sol_complex;
  TrilinosWrappers::MPI::Vector system_rhs_complex;

  sol_complex.reinit(complex_cpu_set, mpi_communicator);
  system_rhs_complex.reinit(complex_cpu_set, mpi_communicator);

  for (auto i : this_cpu_set)
    {
      system_rhs_complex(i)                       = system_rhs(i);
      system_rhs_complex(i + this_cpu_set.size()) = system_rhs_imag(i);
    }

  ConstrainedComplexOperator<TrilinosWrappers::MPI::Vector, BEMProblem<dim>> cc(
    *this, constraints, constraints_imag, constr_cpu_set, mpi_communicator);

  cc.distribute_rhs(system_rhs_complex);
  system_rhs_complex.compress(VectorOperation::insert);

  if (solution_method == "Direct")
    {
      // TODO: find a way to assemble a double size preconditioner
      TrilinosWrappers::SparseMatrix    band_system_complex;
      TrilinosWrappers::PreconditionILU preconditioner_complex;
      TrilinosWrappers::SparsityPattern preconditioner_sparsity_pattern_complex;
      preconditioner_sparsity_pattern_complex.reinit(complex_cpu_set,
                                                     mpi_communicator,
                                                     (types::global_dof_index)
                                                       preconditioner_band);

      for (auto i : this_cpu_set)
        {
          auto from = (i >= preconditioner_band / 2) ?
                        (i - preconditioner_band / 2) :
                        (types::global_dof_index)0;
          auto to =
            std::min((types::global_dof_index)(i + preconditioner_band / 2),
                     this_cpu_set.size());

          for (auto j = from; j < to; ++j)
            {
              preconditioner_sparsity_pattern_complex.add(i, j);
              preconditioner_sparsity_pattern_complex.add(
                i + this_cpu_set.size(), j + this_cpu_set.size());
            }
        }
      preconditioner_sparsity_pattern_complex.compress();
      band_system_complex.reinit(preconditioner_sparsity_pattern_complex);

      for (auto i : this_cpu_set)
        {
          if (constraints.is_constrained(i))
            {
              band_system_complex.add(i, i, 1);
              band_system_complex.add(i + this_cpu_set.size(),
                                      i + this_cpu_set.size(),
                                      1);
            }
          else
            {
              auto from = (i >= preconditioner_band / 2) ?
                            (i - preconditioner_band / 2) :
                            (types::global_dof_index)0;
              auto to =
                std::min((types::global_dof_index)i + preconditioner_band / 2,
                         this_cpu_set.size());

              for (auto j = from; j < to; ++j)
                {
                  if (dirichlet_nodes(i) == 0)
                    {
                      band_system_complex.add(i, j, neumann_matrix(i, j));
                      band_system_complex.add(i + this_cpu_set.size(),
                                              j + this_cpu_set.size(),
                                              neumann_matrix(i, j));
                      // with Robin nodes, the only paired value appears at
                      // this_cpu_set.size() columns on the right - this is only
                      // engaged in very small problems
                      // if ((robin_nodes(i) == 1) &&
                      //     (i + this_cpu_set.size() == j))
                      //   {
                      //     band_system_complex.add(i,
                      //                             j,
                      //                             -dirichlet_matrix(i, i) *
                      //                               robin_matrix_diagonal_imag(
                      //                                 i));
                      //     band_system_complex.add(j,
                      //                             i,
                      //                             dirichlet_matrix(i, i) *
                      //                               robin_matrix_diagonal(i));
                      //   }

                      if (i == j)
                        {
                          band_system_complex.add(i, j, alpha(i));
                          band_system_complex.add(i + this_cpu_set.size(),
                                                  j + this_cpu_set.size(),
                                                  alpha(i));

                          if (robin_nodes(i) == 1)
                            {
                              band_system_complex.add(i,
                                                      j,
                                                      dirichlet_matrix(i, j) *
                                                        robin_matrix_diagonal(
                                                          i));
                              band_system_complex.add(
                                i,
                                j,
                                dirichlet_matrix(i, j) *
                                  robin_matrix_diagonal_imag(i));
                            }
                        }
                    }
                  else
                    {
                      band_system_complex.add(i, j, -dirichlet_matrix(i, j));
                      band_system_complex.add(i + this_cpu_set.size(),
                                              j + this_cpu_set.size(),
                                              -dirichlet_matrix(i, j));
                    }
                }
            }
        }
      preconditioner_complex.initialize(band_system_complex);

      sol_complex.sadd(1., 0., system_rhs_complex);
      solver.solve(cc, sol_complex, system_rhs_complex, preconditioner_complex);
    }
  else
    {
      AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));
      AssertThrow(false, ExcMessage("Unimplemented"));

      TrilinosWrappers::PreconditionILU &fma_preconditioner =
        fma.FMA_preconditioner(alpha, constraints);
      solver.solve(cc, sol, system_rhs, fma_preconditioner);
    }

  for (auto i : this_cpu_set)
    {
      if (neumann_nodes(i) == 1)
        {
          phi(i)      = sol_complex(i);
          phi_imag(i) = sol_complex(i + this_cpu_set.size());
        }
      else if (dirichlet_nodes(i) == 1)
        {
          dphi_dn(i)      = sol_complex(i);
          dphi_dn_imag(i) = sol_complex(i + this_cpu_set.size());
        }
      else
        {
          Assert(robin_nodes(i) == 1, ExcInternalError());
          phi(i)      = sol_complex(i);
          phi_imag(i) = sol_complex(i + this_cpu_set.size());
          // retrieval of dphi_dn using complex coefficients
          std::complex<double> c0_c1(robin_matrix_diagonal(i),
                                     robin_matrix_diagonal_imag(i));
          std::complex<double> c2_c1(robin_rhs(i), robin_rhs_imag(i));
          std::complex<double> ph(phi(i), phi_imag(i));
          std::complex<double> dph_dn(c2_c1 - c0_c1 * ph);
          dphi_dn(i)      = std::real(dph_dn);
          dphi_dn_imag(i) = std::imag(dph_dn);
        }
    }

  phi(this_cpu_set.nth_index_in_set(0)) = phi(this_cpu_set.nth_index_in_set(0));
  phi_imag(this_cpu_set.nth_index_in_set(0)) =
    phi_imag(this_cpu_set.nth_index_in_set(0));
  dphi_dn(this_cpu_set.nth_index_in_set(0)) =
    dphi_dn(this_cpu_set.nth_index_in_set(0));
  dphi_dn_imag(this_cpu_set.nth_index_in_set(0)) =
    dphi_dn_imag(this_cpu_set.nth_index_in_set(0));
  phi.compress(VectorOperation::insert);
  phi_imag.compress(VectorOperation::insert);
  dphi_dn.compress(VectorOperation::insert);
  dphi_dn_imag.compress(VectorOperation::insert);
}

// This method performs a Bem resolution,
// either in a direct or multipole method
template <int dim>
void
BEMProblem<dim>::solve(TrilinosWrappers::MPI::Vector &      phi,
                       TrilinosWrappers::MPI::Vector &      dphi_dn,
                       const TrilinosWrappers::MPI::Vector &tmp_rhs,
                       bool                                 reset_matrix)
{
  if (reset_matrix)
    {
      if (solution_method == "Direct")
        {
          assemble_system();
        }
      else
        {
          AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));

          fma.generate_octree_blocking();
          fma.direct_integrals();
          fma.multipole_integrals();
        }
    }

  solve_system(phi, dphi_dn, tmp_rhs);
}

template <int dim>
void
BEMProblem<dim>::solve(TrilinosWrappers::MPI::Vector &      phi,
                       TrilinosWrappers::MPI::Vector &      phi_imag,
                       TrilinosWrappers::MPI::Vector &      dphi_dn,
                       TrilinosWrappers::MPI::Vector &      dphi_dn_imag,
                       const TrilinosWrappers::MPI::Vector &tmp_rhs,
                       const TrilinosWrappers::MPI::Vector &tmp_rhs_imag,
                       bool                                 reset_matrix)
{
  if (reset_matrix)
    {
      if (solution_method == "Direct")
        {
          assemble_system();
        }
      else
        {
          AssertThrow(dim == 3, ExcMessage("FMA only works in 3D"));

          fma.generate_octree_blocking();
          fma.direct_integrals();
          fma.multipole_integrals();
        }
    }

  solve_system(phi, phi_imag, dphi_dn, dphi_dn_imag, tmp_rhs, tmp_rhs_imag);
}

template <int dim>
void
BEMProblem<dim>::compute_constraints(
  IndexSet &                           c_cpu_set,
  AffineConstraints<double> &          c,
  const TrilinosWrappers::MPI::Vector &tmp_rhs)
{
  Teuchos::TimeMonitor LocalTimer(*ConstraintsTime);
  // We need both the normal vector and surface gradients to apply correctly
  // dirichlet-dirichlet double node constraints. compute_normals();
  compute_surface_gradients(tmp_rhs);

  // communication is needed here: there is one matrix per process: thus the
  // vector needed to set inhomogeneities has to be copied locally
  Vector<double> localized_surface_gradients(
    get_vector_surface_gradients_solution());
  Vector<double> localized_normals(vector_normals_solution);
  Vector<double> localized_dirichlet_nodes(dirichlet_nodes);
  Vector<double> localized_neumann_nodes(neumann_nodes);
  Vector<double> localized_robin_nodes(robin_nodes);
  // Vector<double> localized_robin_flags(robin_flags);
  // Vector<double> localized_robin_rhs(robin_rhs);
  // Vector<double> localized_robin_matrix_diagonal(robin_matrix_diagonal);
  Vector<double> loc_tmp_rhs(tmp_rhs.size());
  loc_tmp_rhs = tmp_rhs;

  // we start clearing the constraint matrix
  c.clear();

  // here we prepare the constraint matrix so as to account for the presence
  // hanging nodes
  AffineConstraints<double> c_hn;
  DoFTools::make_hanging_node_constraints(dh, c_hn);
  c_hn.close();

  std::vector<types::subdomain_id> dofs_domain_association(dh.n_dofs());
  DoFTools::get_subdomain_association(dh, dofs_domain_association);
  // here we prepare the constraint matrix so as to account for the presence of
  // double and triple dofs

  // we start looping on the dofs
  for (types::global_dof_index i = 0; i < tmp_rhs.size(); i++)
    {
      // in the next line we compute the "first" among the set of double nodes:
      // this node is the first dirichlet node in the set, and if no dirichlet
      // node is there, we get the first neumann node
      auto doubles        = double_nodes_set[i];
      auto firstOfDoubles = *doubles.begin();
      for (auto j : doubles)
        {
          // if(this_cpu_set.is_element(j))
          if (localized_dirichlet_nodes(j) == 1)
            {
              firstOfDoubles = j;
              break;
            }
        }
      // do not bind from a robin node if possible
      if (localized_robin_nodes(firstOfDoubles) == 1)
        {
          // only neumann and robins remain
          for (auto j : doubles)
            {
              // if(this_cpu_set.is_element(j))
              if (localized_neumann_nodes(j) == 1)
                {
                  firstOfDoubles = j;
                  break;
                }
            }
        }

      // for each set of double nodes, we will perform the correction only once,
      // and precisely when the current node is the first of the set
      if (i == firstOfDoubles)
        {
          // std::string type = "dirichlet";
          // if (localized_neumann_nodes(i))
          //   {
          //     type = "neumann";
          //   }
          // if (localized_robin_nodes(i))
          //   {
          //     type = "robin";
          //   }
          // pcout << "processing constraints from " << i << " of type " << type
          //       << " coincident with other " << (doubles.size() - 1)
          //       << std::endl;

          // if (!localized_robin_nodes(i) && localized_robin_flags(i))
          //   {
          //     if (!c.is_constrained(i))
          //       {
          //         c.add_line(i);
          //         if (localized_dirichlet_nodes(i))
          //           {
          //             c.set_inhomogeneity(
          //               i,
          //               localized_robin_rhs(i) -
          //                 loc_tmp_rhs(i) *
          //                 localized_robin_matrix_diagonal(i));
          //           }
          //         else
          //           {
          //             c.set_inhomogeneity(i,
          //                                 (localized_robin_rhs(i) -
          //                                  loc_tmp_rhs(i)) /
          //                                   localized_robin_matrix_diagonal(i));
          //           }
          //         pcout
          //           << "node " << i << " of type " << type
          //           << " is on the boundary with a robin condition -> set
          //           inhomogeneity "
          //           << c.get_inhomogeneity(i)
          //           << " imposed value from own condition is " <<
          //           loc_tmp_rhs(i)
          //           << std::endl;
          //       }
          //   }

          // the vector entry corresponding to the first node of the set does
          // i is the source of the constraints, thus we erase it from the set
          doubles.erase(i);

          // TODO: when coinciding with Robin nodes, the rhs should be updated
          // to reflect the known value
          // will these global constraints solve the need for then
          // redistributing the updated rhs?

          // if the current (first) node is a dirichlet node, for all its
          // neumann doubles we will impose that the potential is equal to that
          // of the first node: this means that in the matrix vector product we
          // will put the potential value of the double node
          if (localized_dirichlet_nodes(i) == 1)
            {
              for (auto j : doubles)
                {
                  // if(this_cpu_set.is_element(j))
                  {
                    if (localized_dirichlet_nodes(j) == 1)
                      {
                        // this is the dirichlet-dirichlet case on flat edges:
                        // here we impose that dphi_dn on the two (or more)
                        // sides is equal.
                        double normal_distance = 0;

                        // types::global_dof_index owner_el_1 =
                        // DoFTools::count_dofs_with_subdomain_association (dh,
                        // dofs_domain_association[i]); types::global_dof_index
                        // owner_el_2 =
                        // DoFTools::count_dofs_with_subdomain_association (dh,
                        // dofs_domain_association[*it]);

                        for (unsigned int idim = 0; idim < dim; ++idim)
                          {
                            types::global_dof_index dummy_1 =
                              sub_wise_to_original[i];
                            types::global_dof_index dummy_2 =
                              sub_wise_to_original[j];

                            types::global_dof_index index1 =
                              vec_original_to_sub_wise[gradient_dh.n_dofs() /
                                                         dim * idim +
                                                       dummy_1];
                            types::global_dof_index index2 =
                              vec_original_to_sub_wise[gradient_dh.n_dofs() /
                                                         dim * idim +
                                                       dummy_2];

                            normal_distance += localized_normals[index1] *
                                               localized_normals[index2];
                          }

                        // TODO: validate
                        normal_distance /= normal_distance;
                        if (normal_distance < 1e-4)
                          {
                            c.add_line(j);
                            c.add_entry(j, i, 1);
                          }
                        else if (continuos_gradient)
                          {
                            // this is the dirichlet-dirichlet case on sharp
                            // edges: both normal gradients can be computed from
                            // surface gradients of phi and assingned as BC
                            double norm_i_norm_j = 0;
                            double surf_j_norm_i = 0;
                            double surf_i_norm_j = 0;

                            // types::global_dof_index owner_el_1 =
                            // DoFTools::count_dofs_with_subdomain_association
                            // (dh, dofs_domain_association[i]);
                            // types::global_dof_index owner_el_2 =
                            // DoFTools::count_dofs_with_subdomain_association
                            // (dh, dofs_domain_association[*it]);

                            // We no longer have a std::vector of Point<dim> so
                            // we need to perform the scalar product
                            for (unsigned int idim = 0; idim < dim; ++idim)
                              {
                                types::global_dof_index dummy_1 =
                                  sub_wise_to_original[i];
                                types::global_dof_index dummy_2 =
                                  sub_wise_to_original[j];

                                types::global_dof_index index1 =
                                  vec_original_to_sub_wise
                                    [gradient_dh.n_dofs() / dim * idim +
                                     dummy_1];
                                types::global_dof_index index2 =
                                  vec_original_to_sub_wise
                                    [gradient_dh.n_dofs() / dim * idim +
                                     dummy_2];

                                norm_i_norm_j += localized_normals[index1] *
                                                 localized_normals[index2];
                                surf_j_norm_i +=
                                  localized_surface_gradients[index2] *
                                  localized_normals[index1];
                                surf_i_norm_j +=
                                  localized_surface_gradients[index1] *
                                  localized_normals[index2];
                              }

                            double this_normal_gradient =
                              (1.0 / (1.0 - pow(norm_i_norm_j, 2))) *
                              (surf_j_norm_i +
                               (surf_i_norm_j) * (norm_i_norm_j));
                            double other_normal_gradient =
                              (1.0 / (1.0 - pow(norm_i_norm_j, 2))) *
                              (surf_i_norm_j +
                               (surf_j_norm_i) * (norm_i_norm_j));

                            c.add_line(i);
                            c.set_inhomogeneity(i, this_normal_gradient);
                            c.add_line(j);
                            c.set_inhomogeneity(j, other_normal_gradient);
                          }
                      }
                    else
                      {
                        c.add_line(j);
                        c.set_inhomogeneity(j, loc_tmp_rhs(i));
                        // pcout << "setting inhomogeneity on robin node " << j
                        //       << " = " << c.get_inhomogeneity(j) <<
                        //       std::endl;

                        // if (!c.is_constrained(i))
                        //   {
                        //     // if j is robin, it could also set
                        //     inhomogenerity on i c.add_line(i);
                        //     c.set_inhomogeneity(
                        //       i,
                        //       (localized_robin_rhs(j) - loc_tmp_rhs(i)) /
                        //         localized_robin_matrix_diagonal(j));
                        //     // pcout << "setting inhomogeneity on dirichlet
                        //     node
                        //     // "
                        //     //       << i << " = " << c.get_inhomogeneity(i)
                        //     //       << std::endl;
                        //   }
                      }
                  }
                }
            }

          // if the current (first) node is a neumann node, for all its doubles
          // we will impose that the potential is equal to that of the first
          // node: this means that in the matrix vector product we will put the
          // difference between the potential at the fist node in the doubles
          // set, and the current double node
          if (localized_dirichlet_nodes(i) == 0)
            {
              for (auto j : doubles)
                {
                  c.add_line(j);
                  c.add_entry(j, i, 1);

                  // // if there's a robin, we can set the inhomogeneity
                  // if (localized_neumann_nodes(i) !=
                  // localized_neumann_nodes(j))
                  //   {
                  //     pcout << "setting inhomogeneity on robin node " << j
                  //           << " = " << c.get_inhomogeneity(j) << std::endl;

                  //     // j is a robin node, can set the rhs from the neumann
                  //     // should in some way consider the normals
                  //     c.add_line(i);
                  //     c.set_inhomogeneity(i,
                  //                         localized_robin_rhs(j) -
                  //                           localized_robin_matrix_diagonal(j)
                  //                           *
                  //                             loc_tmp_rhs(i));
                  //     pcout << "setting inhomogeneity on neumann node " << i
                  //           << " = " << c.get_inhomogeneity(i) << std::endl;
                  //   }
                }
            }
        }

      // pcout << "processed double node constraints for dof " << i <<
      // std::endl;
    }

  c.merge(c_hn);
  c.close();

  c_cpu_set.clear();
  c_cpu_set.set_size(this_cpu_set.size());
  for (auto i : this_cpu_set)
    {
      c_cpu_set.add_index(i);
      if (c.is_constrained(i))
        {
          const std::vector<std::pair<types::global_dof_index, double>>
            *entries = c.get_constraint_entries(i);
          for (const auto &pair : *entries)
            {
              c_cpu_set.add_index(pair.first);
            }
        }
    }
  c_cpu_set.compress();
}

template <int dim>
void
BEMProblem<dim>::assemble_preconditioner()
{
  if (!is_preconditioner_initialized)
    {
      for (auto i : this_cpu_set)
        {
          types::global_dof_index start_helper =
            ((i) > preconditioner_band / 2) ? (i - preconditioner_band / 2) :
                                              ((types::global_dof_index)0);
          for (types::global_dof_index j = start_helper;
               j <
               std::min((types::global_dof_index)(i + preconditioner_band / 2),
                        (types::global_dof_index)dh.n_dofs());
               ++j)
            {
              preconditioner_sparsity_pattern.add(i, j);
            }
        }
      preconditioner_sparsity_pattern.compress();
      band_system.reinit(preconditioner_sparsity_pattern);
      is_preconditioner_initialized = true;
    }
  else
    {
      band_system = 0;
    }

  for (auto i : this_cpu_set)
    {
      if (constraints.is_constrained(i))
        {
          band_system.add(i, i, 1);
        }

      types::global_dof_index start_helper = ((i) > preconditioner_band / 2) ?
                                               (i - preconditioner_band / 2) :
                                               ((types::global_dof_index)0);

      for (types::global_dof_index j = start_helper;
           j < std::min((types::global_dof_index)i + preconditioner_band / 2,
                        (types::global_dof_index)dh.n_dofs());
           ++j)
        {
          if (!constraints.is_constrained(i))
            {
              if (dirichlet_nodes(i) == 0)
                {
                  // Nodo di Neumann - or Robin
                  band_system.add(i, j, neumann_matrix(i, j));

                  if (i == j)
                    {
                      band_system.add(i, j, alpha(i));
                      // TODO: account for Robin node
                      if (robin_nodes(i) == 1)
                        {
                          band_system.add(
                            i, j, dirichlet_matrix(i, j) * robin_rhs(i));
                        }
                    }
                }
              else
                {
                  // Nodo di Dirichlet
                  band_system.add(i, j, -dirichlet_matrix(i, j));
                }
            }
        }
    }

  preconditioner.initialize(band_system);
}

template <int dim>
void
BEMProblem<dim>::compute_gradients(
  const TrilinosWrappers::MPI::Vector &glob_phi,
  const TrilinosWrappers::MPI::Vector &glob_dphi_dn)
{
  Teuchos::TimeMonitor LocalTimer(*GradientTime);

  // We need the solution to be stored on a parallel vector with ghost
  // elements. We let Trilinos take care of it.
  TrilinosWrappers::MPI::Vector phi(ghosted_set);
  phi.reinit(glob_phi, false, true);
  TrilinosWrappers::MPI::Vector dphi_dn(ghosted_set);
  dphi_dn.reinit(glob_dphi_dn, false, true);

  // We reinit the gradient solution
  get_vector_gradients_solution().reinit(vector_this_cpu_set, mpi_communicator);

  typedef typename DoFHandler<dim - 1, dim>::active_cell_iterator cell_it;

  // The matrix and rhs of our problem. We must decide if compute the mass
  // matrix just once and for all or not.
  TrilinosWrappers::SparseMatrix vector_gradients_matrix;
  TrilinosWrappers::MPI::Vector  vector_gradients_rhs(vector_this_cpu_set,
                                                     mpi_communicator);
  vector_gradients_matrix.reinit(vector_sparsity_pattern);

  // The vector FEValues to used in the assemblage
  FEValues<dim - 1, dim> vector_fe_v(*mapping,
                                     *gradient_fe,
                                     *quadrature,
                                     update_values | update_gradients |
                                       update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

  // The scalar FEValues to interpolate the known value of phi
  FEValues<dim - 1, dim> fe_v(*mapping,
                              *fe,
                              *quadrature,
                              update_values | update_gradients |
                                update_normal_vectors |
                                update_quadrature_points | update_JxW_values);

  const unsigned int vector_n_q_points    = vector_fe_v.n_quadrature_points;
  const unsigned int vector_dofs_per_cell = gradient_fe->dofs_per_cell;
  std::vector<types::global_dof_index> vector_local_dof_indices(
    vector_dofs_per_cell);


  std::vector<Tensor<1, dim>> phi_surf_grads(vector_n_q_points);
  std::vector<double>         phi_norm_grads(vector_n_q_points);
  std::vector<Vector<double>> q_vector_normals_solution(vector_n_q_points,
                                                        Vector<double>(dim));

  FullMatrix<double> local_gradients_matrix(vector_dofs_per_cell,
                                            vector_dofs_per_cell);
  Vector<double>     local_gradients_rhs(vector_dofs_per_cell);

  std::vector<Point<dim>> support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*mapping,
                                                     dh,
                                                     support_points);
  std::vector<types::global_dof_index> face_dofs(fe->dofs_per_face);

  Quadrature<dim - 1>    dummy_quadrature(fe->get_unit_support_points());
  FEValues<dim - 1, dim> dummy_fe_v(*mapping,
                                    *fe,
                                    dummy_quadrature,
                                    update_values | update_gradients |
                                      update_normal_vectors |
                                      update_quadrature_points);

  const unsigned int                   dofs_per_cell = fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int          n_q_points = dummy_fe_v.n_quadrature_points;
  std::vector<Tensor<1, dim>> dummy_phi_surf_grads(n_q_points);

  cell_it vector_cell = gradient_dh.begin_active();
  cell_it cell = dh.begin_active(), endc = dh.end();
  for (; cell != endc; ++cell, ++vector_cell)
    {
      Assert(cell->index() == vector_cell->index(), ExcInternalError());
      Assert(cell->subdomain_id() == vector_cell->subdomain_id(),
             ExcInternalError());

      if (cell->subdomain_id() == this_mpi_process)
        {
          fe_v.reinit(cell);
          vector_fe_v.reinit(vector_cell);

          local_gradients_matrix = 0;
          local_gradients_rhs    = 0;

          const std::vector<Tensor<1, dim>> &vector_node_normals =
            vector_fe_v.get_normal_vectors();
          fe_v.get_function_gradients(phi, phi_surf_grads);
          fe_v.get_function_values(dphi_dn, phi_norm_grads);
          unsigned int comp_i, comp_j;

          for (unsigned int q = 0; q < vector_n_q_points; ++q)
            {
              Tensor<1, dim> node_normal_grad_dir;
              for (unsigned int i = 0; i < dim; ++i)
                {
                  node_normal_grad_dir[i] = q_vector_normals_solution[q][i];
                }
              Tensor<1, dim> gradient =
                vector_node_normals[q] * phi_norm_grads[q] + phi_surf_grads[q];
              for (unsigned int i = 0; i < vector_dofs_per_cell; ++i)
                {
                  comp_i = gradient_fe->system_to_component_index(i).first;
                  for (unsigned int j = 0; j < vector_dofs_per_cell; ++j)
                    {
                      comp_j = gradient_fe->system_to_component_index(j).first;
                      if (comp_i == comp_j)
                        {
                          local_gradients_matrix(i, j) +=
                            vector_fe_v.shape_value(i, q) *
                            vector_fe_v.shape_value(j, q) * vector_fe_v.JxW(q);
                        }
                    }

                  local_gradients_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                            gradient[comp_i] *
                                            vector_fe_v.JxW(q);
                }
            }
          vector_cell->get_dof_indices(vector_local_dof_indices);

          vector_constraints.distribute_local_to_global(
            local_gradients_matrix,
            local_gradients_rhs,
            vector_local_dof_indices,
            vector_gradients_matrix,
            vector_gradients_rhs);
        }
    }

  // At this point we can compress everything and solve via GMRES.
  vector_gradients_matrix.compress(VectorOperation::add);
  vector_gradients_rhs.compress(VectorOperation::add);

  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(
    solver_control,
    SolverGMRES<TrilinosWrappers::MPI::Vector>::AdditionalData(1000));

  TrilinosWrappers::PreconditionAMG mass_prec;
  mass_prec.initialize(vector_gradients_matrix);
  solver.solve(vector_gradients_matrix,
               get_vector_gradients_solution(),
               vector_gradients_rhs,
               mass_prec);

  vector_constraints.distribute(get_vector_gradients_solution());
}

template <int dim>
void
BEMProblem<dim>::compute_surface_gradients(
  const TrilinosWrappers::MPI::Vector &tmp_rhs)
{
  Teuchos::TimeMonitor          LocalTimer(*SurfaceGradientTime);
  TrilinosWrappers::MPI::Vector phi(ghosted_set);
  phi.reinit(tmp_rhs, false, true);

  get_vector_surface_gradients_solution().reinit(vector_this_cpu_set,
                                                 mpi_communicator);

  typedef typename DoFHandler<dim - 1, dim>::active_cell_iterator cell_it;

  TrilinosWrappers::SparseMatrix vector_surface_gradients_matrix;
  TrilinosWrappers::MPI::Vector  vector_surface_gradients_rhs(
    vector_this_cpu_set, mpi_communicator);

  vector_surface_gradients_matrix.reinit(vector_sparsity_pattern);

  FEValues<dim - 1, dim> vector_fe_v(*mapping,
                                     *gradient_fe,
                                     *quadrature,
                                     update_values | update_gradients |
                                       update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

  FEValues<dim - 1, dim> fe_v(*mapping,
                              *fe,
                              *quadrature,
                              update_values | update_gradients |
                                update_normal_vectors |
                                update_quadrature_points | update_JxW_values);

  const unsigned int vector_n_q_points    = vector_fe_v.n_quadrature_points;
  const unsigned int vector_dofs_per_cell = gradient_fe->dofs_per_cell;
  std::vector<types::global_dof_index> vector_local_dof_indices(
    vector_dofs_per_cell);

  std::vector<Tensor<1, dim>> phi_surf_grads(vector_n_q_points);
  std::vector<double>         phi_norm_grads(vector_n_q_points);
  std::vector<Vector<double>> q_vector_normals_solution(vector_n_q_points,
                                                        Vector<double>(dim));

  FullMatrix<double> local_gradients_matrix(vector_dofs_per_cell,
                                            vector_dofs_per_cell);
  Vector<double>     local_gradients_rhs(vector_dofs_per_cell);

  std::vector<Point<dim>> support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*mapping,
                                                     dh,
                                                     support_points);
  std::vector<types::global_dof_index> face_dofs(fe->dofs_per_face);

  Quadrature<dim - 1>    dummy_quadrature(fe->get_unit_support_points());
  FEValues<dim - 1, dim> dummy_fe_v(*mapping,
                                    *fe,
                                    dummy_quadrature,
                                    update_values | update_gradients |
                                      update_normal_vectors |
                                      update_quadrature_points);

  const unsigned int                   dofs_per_cell = fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int          n_q_points = dummy_fe_v.n_quadrature_points;
  std::vector<Tensor<1, dim>> dummy_phi_surf_grads(n_q_points);

  cell_it vector_cell = gradient_dh.begin_active();
  cell_it cell = dh.begin_active(), endc = dh.end();
  for (; cell != endc; ++cell, ++vector_cell)
    {
      Assert(cell->index() == vector_cell->index(), ExcInternalError());
      Assert(cell->subdomain_id() == vector_cell->subdomain_id(),
             ExcInternalError());

      if (cell->subdomain_id() == this_mpi_process)
        {
          fe_v.reinit(cell);
          vector_fe_v.reinit(vector_cell);
          local_gradients_matrix = 0;
          local_gradients_rhs    = 0;
          fe_v.get_function_gradients(phi, phi_surf_grads);
          unsigned int comp_i, comp_j;

          for (unsigned int q = 0; q < vector_n_q_points; ++q)
            {
              Tensor<1, dim> gradient = phi_surf_grads[q];
              for (unsigned int i = 0; i < vector_dofs_per_cell; ++i)
                {
                  comp_i = gradient_fe->system_to_component_index(i).first;
                  for (unsigned int j = 0; j < vector_dofs_per_cell; ++j)
                    {
                      comp_j = gradient_fe->system_to_component_index(j).first;
                      if (comp_i == comp_j)
                        {
                          local_gradients_matrix(i, j) +=
                            vector_fe_v.shape_value(i, q) *
                            vector_fe_v.shape_value(j, q) * vector_fe_v.JxW(q);
                        }
                    }

                  local_gradients_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                            gradient[comp_i] *
                                            vector_fe_v.JxW(q);
                }
            }
          vector_cell->get_dof_indices(vector_local_dof_indices);

          vector_constraints.distribute_local_to_global(
            local_gradients_matrix,
            local_gradients_rhs,
            vector_local_dof_indices,
            vector_surface_gradients_matrix,
            vector_surface_gradients_rhs);
        }
    }

  vector_surface_gradients_matrix.compress(VectorOperation::add);
  vector_surface_gradients_rhs.compress(VectorOperation::add);

  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(
    solver_control,
    SolverGMRES<TrilinosWrappers::MPI::Vector>::AdditionalData(1000));

  TrilinosWrappers::PreconditionAMG mass_prec;
  mass_prec.initialize(vector_surface_gradients_matrix);

  solver.solve(vector_surface_gradients_matrix,
               get_vector_surface_gradients_solution(),
               vector_surface_gradients_rhs,
               mass_prec);

  vector_constraints.distribute(get_vector_surface_gradients_solution());
}

template <int dim>
void
BEMProblem<dim>::compute_normals()
{
  Teuchos::TimeMonitor LocalTimer(*NormalsTime);
  vector_normals_solution.reinit(vector_this_cpu_set, mpi_communicator);

  typedef typename DoFHandler<dim - 1, dim>::active_cell_iterator cell_it;

  TrilinosWrappers::SparseMatrix vector_normals_matrix;
  TrilinosWrappers::MPI::Vector  vector_normals_rhs(vector_this_cpu_set,
                                                   mpi_communicator);

  vector_normals_matrix.reinit(vector_sparsity_pattern);

  FEValues<dim - 1, dim> vector_fe_v(*mapping,
                                     *gradient_fe,
                                     *quadrature,
                                     update_values | update_gradients |
                                       update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

  const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;

  const unsigned int vector_dofs_per_cell = gradient_fe->dofs_per_cell;

  std::vector<types::global_dof_index> vector_local_dof_indices(
    vector_dofs_per_cell);

  std::vector<Vector<double>> q_vector_normals_solution(vector_n_q_points,
                                                        Vector<double>(dim));

  FullMatrix<double> local_normals_matrix(vector_dofs_per_cell,
                                          vector_dofs_per_cell);
  Vector<double>     local_normals_rhs(vector_dofs_per_cell);

  for (const auto &vector_cell : gradient_dh.active_cell_iterators())
    {
      if (vector_cell->subdomain_id() == this_mpi_process)
        {
          vector_fe_v.reinit(vector_cell);
          local_normals_matrix = 0;
          local_normals_rhs    = 0;
          const std::vector<Tensor<1, dim>> &vector_node_normals =
            vector_fe_v.get_normal_vectors();
          unsigned int comp_i, comp_j;

          for (unsigned int q = 0; q < vector_n_q_points; ++q)
            {
              for (unsigned int i = 0; i < vector_dofs_per_cell; ++i)
                {
                  comp_i = gradient_fe->system_to_component_index(i).first;
                  for (unsigned int j = 0; j < vector_dofs_per_cell; ++j)
                    {
                      comp_j = gradient_fe->system_to_component_index(j).first;
                      if (comp_i == comp_j)
                        {
                          local_normals_matrix(i, j) +=
                            vector_fe_v.shape_value(i, q) *
                            vector_fe_v.shape_value(j, q) * vector_fe_v.JxW(q);
                        }
                    }

                  local_normals_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                          vector_node_normals[q][comp_i] *
                                          vector_fe_v.JxW(q);
                }
            }
          vector_cell->get_dof_indices(vector_local_dof_indices);

          vector_constraints.distribute_local_to_global(
            local_normals_matrix,
            local_normals_rhs,
            vector_local_dof_indices,
            vector_normals_matrix,
            vector_normals_rhs);
        }
    }

  vector_normals_matrix.compress(VectorOperation::add);
  vector_normals_rhs.compress(VectorOperation::add);

  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(
    solver_control,
    SolverGMRES<TrilinosWrappers::MPI::Vector>::AdditionalData(1000));
  TrilinosWrappers::PreconditionAMG mass_prec;
  mass_prec.initialize(vector_normals_matrix);

  solver.solve(vector_normals_matrix,
               vector_normals_solution,
               vector_normals_rhs,
               mass_prec);

  vector_constraints.distribute(vector_normals_solution);
}

template <int dim>
void
BEMProblem<dim>::adaptive_refinement(
  const TrilinosWrappers::MPI::Vector &error_vector)
{
  Vector<float>  estimated_error_per_cell(comp_dom.tria.n_active_cells());
  Vector<double> helper(error_vector);

  KellyErrorEstimator<dim - 1, dim>::estimate(
    *mapping, dh, QGauss<dim - 2>(3), {}, helper, estimated_error_per_cell);

  pgr.mark_cells(estimated_error_per_cell, comp_dom.tria);
  comp_dom.tria.prepare_coarsening_and_refinement();
  comp_dom.tria.execute_coarsening_and_refinement();
}

template class BEMProblem<2>;
template class BEMProblem<3>;
