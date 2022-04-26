root_dir=pwd();
# ==============================================================================
# 0 - PACKAGES
# ==============================================================================
using DataFrames;
using LinearAlgebra;
using JuMP;
using Random;
using Gurobi;
using Graphs;
using GraphPlot;
using Plots;
using Cairo;
using Compose;
using MathOptInterface;
using CSV;

# ==============================================================================
# 1 - PARAMETERS
# ==============================================================================

# Paramètre d'affichage
verbose = false

# Chargement des structures (paramètres, instance, grille, etc.)
println("MST_struct.jl ...");include("$(root_dir)/2_functions/MST_struct.jl");

# Grille de taille Nx*Ny
Nx        = 12
Ny        = 7
# Compacité
beta      = 1
# Nombre de conservation features
spec_nb   = 3
# Graine pseudo-aléatoire
rand_seed = 1
# Méthode de résolution
is_beta        = false
is_non_reserve = true
is_callbacks   = true
is_damier      = true
is_rmax        = true
is_decompose   = true

# Divers
if is_rmax
    Rmax = 3
else
    Rmax = Nx + Ny
end

#folder = "24x14_CF3_BLM0.5_fromfiles"
folder = "$(Nx)x$(Ny)_CF$(spec_nb)_BLM$(beta)"
params = Parameters(beta, spec_nb, rand_seed, is_non_reserve, is_callbacks, is_damier, is_rmax, is_decompose, Rmax)

# Chargement des fonctions
println("MST_functions.jl ...");include("$(root_dir)/2_functions/MST_functions.jl");

# Données du graphe de la grille
gridgraph = GridGraph(Nx, Ny)

instance = Instance(params,gridgraph)
#data_dir = "$(root_dir)/1_data/24x14_CF3_BLM0.5_fromfiles"
#pu_fname     = "$(data_dir)/pu.csv"
#cf_fname     = "$(data_dir)/cf.csv"
#coords_fname = "$(data_dir)/coords.csv"
#bound_fname  = "$(data_dir)/bound.csv"
#instance = Instance(pu_fname,cf_fname,bound_fname,beta,[0.2,0.2,0.2],gridgraph)

if verbose
    # Affichage des infos du graphe de la grille
    print_grid_graph(gridgraph)
    # Affichage des infos de l'instance
    print_instance(instance,gridgraph)
    # Affichage les infos du damier
    damier = Damier(gridgraph.Nx, gridgraph.Ny)
    print_damier(damier.Black, damier.White)
end


# ==============================================================================
# 2 - MODEL + SOLVING
# ==============================================================================

# Chargement des modèles
println("MST_models.jl ...");include("$(root_dir)/2_functions/MST_models.jl");

res_dir = "$(root_dir)/3_results/$(folder)";
if !isdir(res_dir)
    mkdir(res_dir);
end

if !is_decompose
    my_model = ReserveSiteSelection_SpatialConstraints(instance, gridgraph, params, is_beta, is_non_reserve, is_callbacks, is_damier, is_rmax)
    t1 = time_ns()
    optimize!(my_model)
    t2 = time_ns()
    computation_time = round((t2-t1)/1e9,digits=2)
    nb_variables     = length(my_model.moi_backend.optimizer.model.variable_info)
    nb_lin_con       = length(my_model.moi_backend.optimizer.model.affine_constraint_info)
    println("Nombre de variables : $(nb_variables)")
    println("Nombre de contraintes : $(nb_lin_con)")
    println("Temps de calcul : $(computation_time)s")
    x_opt,z_opt,u_opt,r_opt    = read_info_graph_reserve(my_model,gridgraph,is_beta)
    Perimetre,Cout,Score,Rayon,Rayon_InReserve = print_info_reserve(x_opt,instance,gridgraph)
    title = "Solution_$(computation_time)s_Score$(Score)_Perimetre$(Perimetre)_Cout$(Cout)_$(is_beta)_$(is_non_reserve)_$(is_callbacks)_$(is_damier)_$(is_rmax)";
    if is_non_reserve
        x_opt,v_opt = read_info_graph_nonreserve(my_model,gridgraph)
        visualisation_reserve_nonreserve_graph(x_opt,u_opt,r_opt,v_opt,"/$(res_dir)/$(title).png",gridgraph)
    else
        visualisation_reserve_graph(x_opt,u_opt,r_opt,"/$(res_dir)/$(title).png",gridgraph)
    end
else
    # 1. Build the gridgraph and the instance of every ball that includes a feasible solution
    dmin = shortest_distances(length(gridgraph.Noeuds),gridgraph.Voisins)
    feasibleCenters = Vector{Int}()
    feasibleBallNodes = Vector{Vector{Int}}()
    feasibleBallGraphs = Vector{Vector{Int}}()
    feasibleRentability = Vector{Float64}()
    feasibleBallInstances = Vector{Instance}()
    println("Targets = $Targets")
    for i in 1:1 #gridgraph.Noeuds
        # create one gridgraph and one instance for each node and set it as the center of the reserve
        ballNodes, ballGraph = get_subgraph_from_center(gridgraph, i, Rmax, dmin)
        # print(ballNodes)
        # print(ballGraph)
        # Quantité de chaque élément i dans chaque noeud j
        N_cf = instance.N_cf
        ballAmount = zeros(N_cf, length(ballNodes))
        for i in 1:N_cf
            ballAmount[i,:] = instance.Amount[i,ballNodes]
        end
        totalAmount = sum(ballAmount, dims=2)
        Targets = instance.Targets

        # test the capacity of the ball to satisfy all constraints
        isFeasible = true
        println("Total amout in ball = $(totalAmount[:,1])")
        for i in 1:N_cf
            if totalAmount[i,1] < Targets[i]
                isFeasible = false
                break
            end
        end

        if isFeasible
            push!(feasibleCenters, i)
            push!(feasibleBallNodes, ballNodes)
            push!(feasibleBallGraphs, ballGraph)
            push!(feasibleRentability, sum(instance.Rentability[ballNodes]))

             # define the remaining data of the subgraph
            ballCost = instance.Cost[ballNodes]
            ballLockedOut = instance.IsLockedOut[ballNodes]
            ballRentability = instance.Rentability[ballNodes]
            ballBoundaryLength = Dict{Pair{Int,Int},Int}()
            for d in ballGraph.Arcs
                ballBoundaryLength[d] = instance.BoundaryLength[ballNodes[d[1]]=>ballNodes[d[2]]]
            end
            ballBoundaryCorrection = zeros(length(ballNodes))
            for i in ballGraph.NoeudsPeripheriques
                ballBoundaryCorrection[i] = 4 - length(ballGraph.Voisins[i])
            end
            ballInstance = Instance(length(ballNodes),instance.N_cf,instance.N_bd,ballGraph.Noeuds,instance.ConservationFeatures,ballCost,ballAmount,Targets,ballRentability,ballBoundaryLength,ballBoundaryCorrection,instance.Beta,ballLockedOut)
            push!(feasibleBallInstances, ballInstance)
        end        
    end

    # 2. Sort the balls by decreasing rentability to start with those that are probably the best 
    sortedFeasible = sortperm(Rentability; rev=true)
    
    # 3. Solve the problem in each ball
    for i in sortedFeasible
        my_model = ReserveSiteSelection_SpatialConstraints(feasibleBallInstances[i], feasibleBallGraphs[i], params, is_beta, is_non_reserve, is_callbacks, is_damier, is_rmax)
        t1 = time_ns()
        optimize!(my_model)
        t2 = time_ns()
        computation_time = round((t2-t1)/1e9,digits=2)
        nb_variables     = length(my_model.moi_backend.optimizer.model.variable_info)
        nb_lin_con       = length(my_model.moi_backend.optimizer.model.affine_constraint_info)
        println("Nombre de variables : $(nb_variables)")
        println("Nombre de contraintes : $(nb_lin_con)")
        println("Temps de calcul : $(computation_time)s")
        x_opt,z_opt,u_opt,r_opt    = read_info_graph_reserve(my_model,gridgraph,is_beta)
    end 
end

