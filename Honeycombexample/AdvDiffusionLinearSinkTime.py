
# import pdb; pdb.set_trace()

# Add Python bindings directory to PATH
import sys, os, scipy, numpy


#ys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

def read_ansys_mesh(mesh_dir, filename, nodes_subset=[], elem_subset=[], debug=False):

    # Load ansys .in file
    f = open(os.path.join(mesh_dir, filename), 'r')
    lines = f.readlines()
    num_lines = len(lines)

    # Initialising empty arrays in which to store node values and coordinates
    nodes_per_elem = 4  # for a tet mesh

    node_array = numpy.empty((0, 1), int)
    node_coordinates = numpy.empty((0, 3), int)
    element_array = numpy.empty((0, 1), int)
    element_nodes_array = numpy.empty((0, nodes_per_elem), int)

    # Add nodes
    for line_idx, line in enumerate(lines):
        if line.split(' ,')[0] == 'NBLOCK':
            for node_line_idx in range(line_idx + 2, num_lines + 1):
                node_line = lines[node_line_idx]
                if node_line.split()[0] == 'N':
                    break
                else:
                    coordinates = node_line.split('       ')[-1]
                    x = float(coordinates[1:17])

                    y = float(coordinates[17:33])
                    z = float(coordinates[33:-1])
                    node_num = int(node_line.split()[0])
                    if node_num in nodes_subset or nodes_subset == []:
                        # Save node numbers (node_num) and coordinates (x, y, z) to arrays
                        node_array = numpy.append(node_array, node_num)
                        node_coordinates = numpy.append(node_coordinates, numpy.array([[x, y, z]]), axis=0)
            break
            print 'X Y Z step'

    # Add elements
    for line_idx, line in enumerate(lines):
        if line.split(' ,')[0] == 'EBLOCK':
            for node_line_idx in range(line_idx + 2, num_lines + 1):
                node_line = lines[node_line_idx]
                if node_line.split() == []:
                    break
                else:
                    element_nodes = node_line.split()[11:-1]
                    element_nodes, idx_array = scipy.unique(scipy.array([int(node) for node in element_nodes]),
                                                            return_index=True)
                    idx_array = [3 if idx == 4 else idx for idx in idx_array]

                    # Reordering the node arrangement
                    renumbered_nodes = scipy.copy(element_nodes)
                    for position, idx in enumerate(idx_array):
                        renumbered_nodes[idx] = element_nodes[position]

                    element_num = int(node_line.split()[10])
                    if element_num in elem_subset or elem_subset == []:
                        # Save element number (element_num) and element nodes (element_nodes) to arrays
                        element_array = numpy.append(element_array, element_num)
                        element_nodes_array = numpy.append(element_nodes_array, numpy.array([renumbered_nodes]), axis=0)
            break
            print 'elem step'

    inlet_node_array = numpy.empty((0, 1), int)
#   outlet_node_array = numpy.empty((0, 1), int)

    # Find which nodes are part of the inlet
    #------------------------------------------
    for line_idx, line in enumerate(lines):
        if line.split(',')[0] == 'CMBLOCK' and line.split(',')[1] == 'INLET':
            for node_line_idx in range(line_idx + 2, num_lines + 1):
                node_line = lines[node_line_idx]
                if node_line.split(',')[0] == 'CMBLOCK':
                    break
                else:
                    inlet_node_row = node_line.split()
                    inlet_node_row = scipy.array([int(node) for node in inlet_node_row])
                    print ' step', inlet_node_row
                    for inlet_node in inlet_node_row:
                        inlet_node_array = numpy.append(inlet_node_array, inlet_node)
            break

    # Find which nodes are part of the outlet
    #------------------------------------------
#   for line_idx, line in enumerate(lines):
#       if line.split(',')[0] == 'CMBLOCK' and line.split(',')[1] == 'OUT2':
#           for node_line_idx in range(line_idx + 2, num_lines + 1):
#               node_line = lines[node_line_idx]
#               if node_line.split(',')[0] == 'CMBLOCK':
#                   break
#               else:
#                   outlet2_node_row = node_line.split()
#                   outlet2_node_row = scipy.array([int(node) for node in outlet2_node_row])
#                   for outlet2_node in outlet1_node_row:
#                       outlet2_node_array = numpy.append(outlet2_node_array, outlet2_node)
#           break

    # Find which nodes are part of the outlet
    #------------------------------------------
#   for line_idx, line in enumerate(lines):
#       if line.split(',')[0] == 'CMBLOCK' and line.split(',')[1] == 'OUT3':
#           for node_line_idx in range(line_idx + 2, num_lines + 1):
#               node_line = lines[node_line_idx]
#               if node_line.split(',')[0] == 'CMBLOCK':
#                   break
#               else:
#                   outlet3_node_row = node_line.split()
#                   outlet3_node_row = scipy.array([int(node) for node in outlet3_node_row])
#                   for outlet3_node in outlet3_node_row:
#                       outlet3_node_array = numpy.append(outlet3_node_array, outlet3_node)
#           break


    # Return node number and coordinate arrays, and element number and element node arrays
    # ------------------------------------------------------------------------------------
#return node_array, node_coordinates, element_array, element_nodes_array, inlet_node_array, outlet2_node_array,outlet3_node_array
    return node_array, node_coordinates, element_array, element_nodes_array, inlet_node_array




# Call the code which reads the ansys mesh
#----------------------------------------
### [node_array, node_coordinates, element_array, element_nodes_array, inlet_node_array,outlet2_node_array,

[node_array, node_coordinates, element_array, element_nodes_array, 
 inlet_node_array] = read_ansys_mesh('./', 'project1Coarse11.in')


# Changing the values in each array to 32 bit integers
node_array = node_array.astype(numpy.int32)
element_array = element_array.astype(numpy.int32)
element_nodes_array = element_nodes_array.astype(numpy.int32)


# Initialise OpenCMISS-Iron
from opencmiss.iron import iron


# Set problem parameters
diff_coeff = 20.  # from readings
initial_conc = 0.001
start_time = 0.0
end_time = 0.1
time_step = 0.001
screen_output_freq = 1  # how many time steps between outputs to screen

# Set problem parameters

(coordinateSystemUserNumber,
 regionUserNumber,
 basisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 geometricFieldUserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 independentFieldUserNumber,
 materialFieldUserNumber,
 equationsSetUserNumber,
 problemUserNumber) = range(1, 14)

#sourceFieldUserNumber,

#ron.DiagnosticsSetOn(iron.DiagnosticTypes.IN, [1, 2, 3, 4, 5], "Diagnostics",
#                     ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Get the computational nodes information
#----------------------------------------
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

number_of_dimensions = 3
number_of_mesh_components = 1
total_number_of_elements = len(element_array)
total_number_of_nodes = len(node_array)
mesh_component_number = 1
nodes_per_elem = 4  # for a tet mesh


# Create a RC coordinate system
#--------------------------------
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()


# Create a region
#------------------------------
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.label = "DiffusionRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear simplex basis
#------------------------------------
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.SIMPLEX)
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX] * 3
basis.CreateFinish()


# Create a tri-quadratic Lagrange basis
#basis = iron.Basis()
#basis.CreateStart(basisUserNumber)
#basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
#basis.numberOfXi = 3
#basis.interpolationXi = [
#       iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * 3
#basis.quadratureNumberOfGaussXi = [3] * 3
#basis.CreateFinish()


# Start the creation of the imported mesh in the region
# -----------------------------------------------
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber, region, number_of_dimensions)
mesh.NumberOfComponentsSet(number_of_mesh_components)
mesh.NumberOfElementsSet(total_number_of_elements)


# Define nodes for the mesh
#-------------------------------
nodes = iron.Nodes()
nodes.CreateStart(region, total_number_of_nodes)

nodes.UserNumbersAllSet(node_array)
nodes.CreateFinish()

elements = iron.MeshElements()
elements.CreateStart(mesh, mesh_component_number, basis)


# Set the nodes pertaining to each element
#----------------------------------------
for idx, elem_num in enumerate(element_array):
    elements.NodesSet(idx + 1, element_nodes_array[idx])

# Refers to elements by their user number as described in the original mesh
elements.UserNumbersAllSet(element_array)
elements.CreateFinish()

mesh.CreateFinish()


# Create a decomposition for the mesh
#------------------------------------
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()


# Create a field for the geometry
#-------------------------------
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()


# Update the geometric field parameters
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

for idx, node_num in enumerate(node_array):
    [x, y, z] = node_coordinates[idx]

    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(node_num), 1, x)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(node_num), 2, y)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(node_num), 3, z)

geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)


# Create standard Diffusion equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()

equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
        iron.EquationsSetTypes.ADVECTION_DIFFUSION_EQUATION,
        iron.EquationsSetSubtypes.NO_SOURCE_ADVECTION_DIFFUSION]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()


#Create source field for equation set
#----------------------------------
#sourceField = iron.Field()
#equationsSet.SourceCreateStart(sourceFieldUserNumber, sourceField)
#sourceField.VariableLabelSet(iron.FieldVariableTypes.U, "Source")
#equationsSet.SourceCreateFinish()

#sourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)


# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Concen")
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()


# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()


## This will change the diffusion coeff
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,diff_coeff)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,diff_coeff)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,diff_coeff)


# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, initial_conc)



# Create independent field
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
independentField.VariableLabelSet(iron.FieldVariableTypes.U, "Velocity")
equationsSet.IndependentCreateFinish()


independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,-5.0)
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,0.0)


#for node_num in enumerate(node_array):
#  independentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 42, 0, int(node_num), int(1), )
#  independentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 42, 0, int(node_num), int(2), )
#  independentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 42, 0, int(node_num), int(3), )



##for vnode in vList:
# independentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(vnode[0]), int(1), vnode[1])
# independentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(vnode[0]), int(2), vnode[2])
# independentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, int(vnode[0]), int(3), vnode[3])



# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()


# Create Adv-Diffusion equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
        		iron.ProblemTypes.ADVECTION_DIFFUSION_EQUATION,
        		iron.ProblemSubtypes.NO_SOURCE_ADVECTION_DIFFUSION]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()


# Create control loops
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.TimesSet(start_time,end_time,time_step)
controlLoop.TimeOutputSet(screen_output_freq)
problem.ControlLoopCreateFinish()



# Create problem solver
#solver = iron.Solver()
#linearSolver = iron.Solver()
#problem.SolversCreateStart()
#problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
#solver.outputType = iron.SolverOutputTypes.PROGRESS
#problem.SolversCreateFinish()

## Create solver equations and add equations set to solver equations
#solver = iron.Solver()
#solverEquations = iron.SolverEquations()
#problem.SolverEquationsCreateStart()
#problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
#solver.SolverEquationsGet(solverEquations)
#solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
#equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
#problem.SolverEquationsCreateFinish()


# Create problem solver
dynamicsolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,dynamicsolver)
dynamicsolver.outputType = iron.SolverOutputTypes.PROGRESS
linearSolver=iron.Solver()
dynamicsolver.DynamicLinearSolverGet(linearSolver)
linearSolver.outputType = iron.SolverOutputTypes.NONE
linearSolver.linearType = iron.LinearSolverTypes.ITERATIVE
linearSolver.LinearIterativeMaximumIterationsSet(1000)
problem.SolversCreateFinish()

## Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()


# Create boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)


# Set maximum concentration (1) for nodes at the inlet
#-----------------------------------------------------
for inlet_node in inlet_node_array:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(inlet_node),1,iron.BoundaryConditionsTypes.FIXED,1.0)

#Set minimum concentration (0) for nodes at the outlet
#for outlet_node in outlet_node_array:
#  boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,int(outlet_node),1,iron.BoundaryConditionsTypes.FIXED,0.0)

solverEquations.BoundaryConditionsCreateFinish()



# Manually generating the field values in one single control loop
current_field_array = numpy.zeros(len(node_array))

print 'Time step:', time_step
print 'Start time:', start_time
print 'End time:', end_time

#### controlLoop.TimesSet(start_time, end_time, time_step)

# Solve the problem using the new time loop
problem.Solve()

# Get the field values
for idx, node_num in enumerate(node_array):
    current_field_array[idx] = dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                               iron.FieldParameterSetTypes.VALUES,
                                                                               1, 1,int(node_num), 1)

print 'The field array for end time', end_time, 'is:\n', current_field_array


# Export results
baseName = "AdvDiffusion"
dataFormat = "PLAIN_TEXT"
fml = iron.FieldMLIO()
fml.OutputCreate(mesh, "", baseName, dataFormat)
fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputWrite("AdvDiffusionExample.xml")
fml.Finalise()


# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("HoneycombAdvDiffVel5p0_0p1", "FORTRAN")
#fields.ElementsExport("BranchEmph3_2p4", "FORTRAN")
fields.Finalise()


iron.Finalise()




