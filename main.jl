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
    Nx        = 2*12
    Ny        = 2*7
    # Compacité
    beta      = -1
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

    # Divers
    Rmax = 6

    #folder = "24x14_CF3_BLM0.5_fromfiles"
    folder = "$(Nx)x$(Ny)_CF$(spec_nb)_BLM$(beta)"


    if is_rmax
        params = Parameters(Nx,Ny,beta,spec_nb,rand_seed,is_non_reserve,is_callbacks,is_damier,is_rmax,Rmax)
    else
        params = Parameters(Nx,Ny,beta,spec_nb,rand_seed,is_non_reserve,is_callbacks,is_damier,is_rmax)
    end

    # Chargement des fonctions
    println("MST_functions.jl ...");include("$(root_dir)/2_functions/MST_functions.jl");


    # Données du graphe de la grille
    gridgraph = GridGraph(params)
    Noeuds              = gridgraph.Noeuds
    Aretes              = gridgraph.Aretes
    Voisins             = gridgraph.Voisins
    Arcs                = gridgraph.Arcs
    NoeudsPeripheriques = gridgraph.NoeudsPeripheriques
    NoeudsFictif        = gridgraph.NoeudsFictif
    AretesFictif        = gridgraph.AretesFictif
    VoisinsFictif       = gridgraph.VoisinsFictif
    ArcsFictif          = gridgraph.ArcsFictif
    alpha               = gridgraph.alpha
    N_noeuds            = length(Noeuds)
    N_noeudsfictif      = length(NoeudsFictif)
    N_aretes            = length(Aretes)

    # Affichage des infos du graphe de la grille
    if verbose
        print_grid_graph(gridgraph)
    end

    instance = InstanceFromParams(params,gridgraph)
    #data_dir = "$(root_dir)/1_data/24x14_CF3_BLM0.5_fromfiles"
    #pu_fname     = "$(data_dir)/pu.csv"
    #cf_fname     = "$(data_dir)/cf.csv"
    #coords_fname = "$(data_dir)/coords.csv"
    #bound_fname  = "$(data_dir)/bound.csv"
    #instance = InstanceFromFiles(pu_fname,cf_fname,bound_fname,beta,[0.2,0.2,0.2],gridgraph)

    ConservationFeatures = instance.ConservationFeatures
    Amount               = instance.Amount
    Cost                 = instance.Cost
    Targets              = instance.Targets
    Rentability          = instance.Rentability
    BoundaryLength       = instance.BoundaryLength
    BoundaryCorrection   = instance.BoundaryCorrection

    # Affichage des infos de l'instance
    if verbose
        print_instance(instance,gridgraph)
    end


    # Infos du damier
    damier = Damier(params)
    NoeudsNoirs  = damier.Black
    NoeudsBlancs = damier.White

    # Affichage les infos du damier
    if verbose
        print_damier(damier)
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
