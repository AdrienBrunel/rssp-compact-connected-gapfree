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

const GRB_ENV = Gurobi.Env()

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
params = Parameters(beta, spec_nb, rand_seed, is_non_reserve, is_callbacks, is_damier, is_beta, is_rmax, is_decompose, Rmax)

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


function solve_instance(instance, gridgraph, params, res_dir)
    if !params.is_decompose
        my_model = ReserveSiteSelection_SpatialConstraints(instance, gridgraph, params)
        t1 = time_ns()
        optimize!(my_model)
        if !has_values(my_model)
            println("The problem is not feasible")
        end
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
        feasibleBallGraphs = Vector{GridGraph}()
        feasibleRentability = Vector{Float64}()
        feasibleBallInstances = Vector{Instance}()
        for i in gridgraph.Noeuds
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

            # test the capacity of the ball to satisfy all constraints
            isFeasible = true
            println("Total amout in ball = $(totalAmount[:,1])")
            for i in 1:N_cf
                if totalAmount[i,1] < instance.Targets[i]
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
                ballInstance = Instance(length(ballNodes), instance.N_cf, instance.N_bd, ballGraph.Noeuds, instance.ConservationFeatures, ballCost, ballAmount, instance.Targets, ballRentability, ballBoundaryLength, ballBoundaryCorrection, instance.Beta,ballLockedOut)
                push!(feasibleBallInstances, ballInstance)
            end        
        end

        # 2. Sort the balls by decreasing rentability to start with those that are probably the best 
        sortedFeasible = sortperm(feasibleRentability; rev=true)
        
        # 3. Solve the problem in each ball
        best_objval = Inf
        best_x = zeros(Int, length(gridgraph.Noeuds))
        best_r = zeros(Int, length(gridgraph.Noeuds))
        best_u = Dict{Pair{Int,Int},Int}()
        for i in sortedFeasible
            mip_model = ReserveSiteSelection_SpatialConstraints(feasibleBallInstances[i], feasibleBallGraphs[i], params, findfirst(feasibleBallNodes[i] .== feasibleCenters[i]), best_objval)

            # set gurobi parameters
            set_optimizer_attribute(mip_model, "OutputFlag",Int(verbose))
            set_optimizer_attribute(mip_model, "TimeLimit", 10.0)
            set_optimizer_attribute(mip_model, "Threads", 1)

            #  call Gurobi
            optimize!(mip_model)

            # get the optimal value
            if !has_values(mip_model)
                println("\n***************************************************")
                println("Center $(feasibleCenters[i]): could not improve best solution")
                println("***************************************************\n")
                continue
            end
            objval =JuMP.objective_value(mip_model)
            println("\n***************************************************")
            println("Center $(feasibleCenters[i]): one new solution with value $objval has been found")
            println("***************************************************\n")

            # update incumbent
            if objval < best_objval
                best_objval = objval
                x_opt,z_opt,u_opt,r_opt = read_info_graph_reserve(mip_model, feasibleBallGraphs[i], params.is_beta)
                Perimetre,Cout,Score,Rayon,Rayon_InReserve = print_info_reserve(x_opt, feasibleBallInstances[i], feasibleBallGraphs[i])
                for j in gridgraph.Noeuds
                    best_x[j] = 0
                    best_r[j] = 0
                end
                for d in gridgraph.Arcs
                    best_u[d] = 0
                end
                for j in 1:length(feasibleBallNodes[i])
                    best_x[feasibleBallNodes[i][j]] = x_opt[j]
                end
                best_r[feasibleCenters[i]] = 1
                for d in feasibleBallGraphs[i].Arcs
                    best_u[feasibleBallNodes[i][d[1]]=>feasibleBallNodes[i][d[2]]] = u_opt[d]
                end
                println("new best solution = $(findall(best_x .== 1))")
            end
        end
        # output the best solution among all balls
        println("\n***************************************************")
        println("Best solution is for center $(findfirst(best_r .== 1)) with value $best_objval")
        println("***************************************************\n")
        Perimetre,Cout,Score,Rayon,Rayon_InReserve = print_info_reserve(best_x,instance,gridgraph)
        title = "Solution_0s_Score$(Score)_Perimetre$(Perimetre)_Cout$(Cout)_$(is_beta)_$(is_non_reserve)_$(is_callbacks)_$(is_damier)_$(is_rmax)"
        visualisation_reserve_graph(best_x,best_u,best_r,"/$(res_dir)/$(title).png",gridgraph)
    end
end

