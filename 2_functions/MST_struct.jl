# ==============================================================================
# 1 - PARAMETERS
# ==============================================================================
struct Parameters
    beta::Float64
    spec_nb::Int64
    rand_seed::Int64
    is_non_reserve::Bool
    is_callbacks::Bool
    is_damier::Bool
    is_rmax::Bool
    Rmax::Int64
end



# ==============================================================================
#  DAMIER
# ==============================================================================
struct Damier
    IsBlack::BitVector
    Black::Vector{Int64}
    White::Vector{Int64}
end

# Fonction générant la liste des noeuds noirs et blancs du damier de la grille
function Damier(Nx::Int, Ny::Int)
    # Initialisation
    IsBlack = falses(Nx*Ny)
    Black   = Vector{Int64}()
    White   = Vector{Int64}()
    
    # Damier
    num_ligne = 0
    for i in 1:(Nx*Ny)
        if mod(i,Nx) == 0
            num_colonne = Nx
        else
            num_colonne = mod(i,Nx)
        end
        
        if mod(i,Nx) == 1
            num_ligne = num_ligne + 1
        end
        
        if mod(num_ligne,2) == mod(num_colonne,2)
            IsBlack[i] = true
            push!(Black,i)
        else
            push!(White,i)
        end
    end
    
    return Damier(IsBlack, Black, White)
end

# ==============================================================================
# 2 - GRID GRAPH
# ==============================================================================
"""
    GridGraph

    Structure describing the discretized marine space as a graph

    # Fields
    * `Nx::Int`: Maximum abscissa of a node of the graph
    * `Ny::Int`: Maximum ordinate of a node of the graph
    * `Noeuds::Vector{Int}`:
    * `Voisins::Vector{Vector{Int}}`:
    * `Arcs::Vector{Pair{Int,Int}}`:
    * `NoeudsPeripheriques::Vector{Int}`:
    * `VoisinsFictif::Vector{Vector{Int}}`:
    * `ArcsFictif::Vector{Pair{Int,Int}}`:
    * `alpha::Int64`:
    * `damier::Damier`:
"""
struct GridGraph
    Nx::Int
    Ny::Int
    Noeuds::Vector{Int}
    Voisins::Vector{Vector{Int}}
    Arcs::Vector{Pair{Int,Int}}
    NoeudsPeripheriques::Vector{Int}
    VoisinsFictif::Vector{Vector{Int}}
    ArcsFictif::Vector{Pair{Int,Int}}
    alpha::Int64
    damier::Damier
end

# Fonction donnant les éléments du graphe de la grille
function GridGraph(Nx::Int, Ny::Int)
    # Noeuds du graphe associé à la grille de taille Nx*Ny
    Noeuds = 1:Nx*Ny;
    
    # Noeud fictif
    alpha = Nx*Ny+1;
    
    # Liste des voisins pour chaque noeud de la grille de taille Nx*Ny
    Voisins = Vector{Vector{Int}}()
    for _ in Noeuds
        push!(Voisins, Vector{Int}())
    end
    for i in Noeuds
        # arêtes horizontales
        if mod(i,Nx) != 0
            push!(Voisins[i], i+1)
            push!(Voisins[i+1], i)
        end
        # arêtes verticales
        if i <= (Ny-1)*Nx
            push!(Voisins[i], i+Nx)
            push!(Voisins[i+Nx], i)
        end
    end
            
    # Arcs du graphe associé à la grille de taille Nx*Ny
    Arcs = Vector{Pair{Int,Int}}()
    for i in Noeuds
        for j in Voisins[i]
            push!(Arcs,i=>j)
        end
    end
    
    # Noeuds périphériques du graphe
    NoeudsPeripheriques = Vector{Int64}()
    for i in Noeuds
        if length(Voisins[i]) < 4
            push!(NoeudsPeripheriques,i)
        end
    end
    
    # Arêtes du graphe associé à la grille de taille Nx*Ny
    ArcsFictif = Vector{Pair{Int,Int}}()
    for i in NoeudsPeripheriques
        # arêtes noeuds fictif avec les noeuds périphériques
        push!(ArcsFictif,(alpha=>i))
        push!(ArcsFictif,(i=>alpha))
    end
    
    # Liste des voisins pour chaque noeud y compris le noeud fictif
    VoisinsFictif = Vector{Vector{Int}}()
    for i in Noeuds
        push!(VoisinsFictif, Vector{Int}())
        for j in  Voisins[i]
            push!(VoisinsFictif[i], j)
        end
    end
    push!(VoisinsFictif, Vector{Int}())
    for i in NoeudsPeripheriques
        push!(VoisinsFictif[i], alpha)
        push!(VoisinsFictif[alpha], i)
    end

    
    # Noeuds noirs et blanc du damier
    damier = Damier(Nx, Ny)
    
    return GridGraph(Nx, Ny, Noeuds,Voisins,Arcs,NoeudsPeripheriques,VoisinsFictif,ArcsFictif,alpha, damier)
end


# ==============================================================================
# 3 - INSTANCE
# ==============================================================================
struct Instance
    N_pu::Int64 #
    N_cf::Int64 #
    N_bd::Int64 #
    PlanningUnits::Vector{Int64}
    ConservationFeatures::Vector{Int64}
    Cost::Vector{Float64}
    Amount::Matrix{Float64}
    Targets::Vector{Float64}
    Rentability::Vector{Float64}
    BoundaryLength::Dict{Pair{Int,Int},Int}
    BoundaryCorrection::Vector{Float64}
    Beta::Float64
    IsLockedOut::Array{Int8,2}
end

# Fonction d'habillage de l'instance
function Instance(params,gridgraph)
    
    # Parameters
    beta      = params.beta
    N_cf      = params.spec_nb
    rand_seed = params.rand_seed
    
    # Données du graphe de la grille
    N_pu                = gridgraph.Nx * gridgraph.Ny
    Noeuds              = gridgraph.Noeuds
    Voisins             = gridgraph.Voisins
    Arcs                = gridgraph.Arcs
    NoeudsPeripheriques = gridgraph.NoeudsPeripheriques
    N_noeuds = length(Noeuds)
    PlanningUnits = Noeuds
    N_bd = length(Arcs)
    
    # Quantité de chaque élément i dans chaque noeud j
    Random.seed!(rand_seed)
    ConservationFeatures = collect(1:N_cf)
    Amount = zeros(N_cf,N_noeuds)
    for i in ConservationFeatures
        for j in Noeuds
            Amount[i,j] = rand(0:5)
        end
    end
    
    # Coût de chaque noeud j
    Cost = zeros(N_noeuds)
    IsLockedOut = zeros(N_pu,1)
    for j in Noeuds
        Cost[j] = rand(1:5)
    end
    
    # La réserve doit recouvrir chaque élément i dans le respect des quantités cibles
    Targets = zeros(N_cf)
    for i in ConservationFeatures
        Targets[i] = 0.2*sum(Amount[i,:])
    end
    
    # Calcul de la rentabilité d'une cellule
    Rentability = zeros(N_noeuds)
    for j in Noeuds
        Rentability[j] = sum(Amount[i,j]/Targets[i] for i in ConservationFeatures)/Cost[j]
    end
    
    # Longueur de chaque arcs
    BoundaryLength = Dict{Pair{Int,Int},Int}()
    for d in Arcs
        BoundaryLength[d] = 1
    end
    
    # Correction
    BoundaryCorrection = zeros(N_noeuds)
    for j in NoeudsPeripheriques
        BoundaryCorrection[j] = 4 - length(Voisins[j])
    end
    
    return Instance(N_pu,N_cf,N_bd,PlanningUnits,ConservationFeatures,Cost,Amount,Targets,Rentability,BoundaryLength,BoundaryCorrection,beta,IsLockedOut)
    
end

function Instance(pu_fname,cf_fname,bound_fname,beta,targets,gridgraph)
        
    # Check if the input data files exist
    for f in [pu_fname,cf_fname,bound_fname]
        if !isfile(f)
            println("WARNING! File $(f) is missing")
        end
    end
    
    # Data of the graph associated with the grid
    Voisins             = gridgraph.Voisins
    NoeudsPeripheriques = gridgraph.NoeudsPeripheriques
    
    # Read input files
    pu_data    = CSV.read(pu_fname, header=1, delim=",")
    cf_data    = CSV.read(cf_fname, header=1, delim=",")
    bound_data = CSV.read(bound_fname, header=1, delim=",")
    
    # Problem size
    N_pu = length(pu_data.id)
    N_cf = length(targets)
    N_bd = length(bound_data.boundary)
    
    # Problem ranges
    PlanningUnits = collect(1:N_pu)
    ConservationFeatures = collect(1:N_cf)
    
    # Cost and status of planning units
    Cost = zeros(N_pu)
    IsLockedOut = zeros(N_pu)
    #LockedOut = Vector{Int64}()
    for j in 1:N_pu
        Cost[j] = pu_data.cost[j]
        IsLockedOut[j] = pu_data.is_locked_out[j]
    end
    
    # Amount and targets of conservation features
    Amount = zeros(N_cf,N_pu)
    Targets = zeros(N_cf)
    for i in 1:N_cf
        Amount[i,1:N_pu] = cf_data[:,i+1]
        Targets[i] = targets[i] * sum(Amount[i,1:N_pu] .* (1 .- IsLockedOut))
    end
    
    # Calcul de la rentabilité d'une cellule
    Rentability = zeros(N_noeuds)
    for j in Noeuds
        Rentability[j] = sum(Amount[i,j]/Targets[i] for i in ConservationFeatures)/Cost[j]
    end
    
    # Boundary length of vertices between two nodes
    BoundaryLength = Dict{Pair{Int,Int},Int}()
    for i in 1:N_bd
        BoundaryLength[bound_data.id1[i]=>bound_data.id2[i]] = bound_data.boundary[i]
        BoundaryLength[bound_data.id2[i]=>bound_data.id1[i]] = bound_data.boundary[i]
    end
    
    # Correction
    BoundaryCorrection = zeros(N_pu)
    for j in NoeudsPeripheriques
        BoundaryCorrection[j] = 4 - length(Voisins[j])
    end
    
    # constructor
    return Instance(N_pu,N_cf,N_bd,PlanningUnits,ConservationFeatures,Cost,Amount,Targets,Rentability,BoundaryLength,BoundaryCorrection,beta,IsLockedOut)
    
end

