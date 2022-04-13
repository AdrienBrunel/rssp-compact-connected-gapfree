# ==============================================================================
# 1 - PARAMETERS
# ==============================================================================
    struct Parameters
        Nx::Int64
        Ny::Int64
        beta::Float64
        spec_nb::Int64
        rand_seed::Int64
        is_non_reserve::Bool
        is_callbacks::Bool
        is_damier::Bool
        is_rmax::Bool
        Rmax::Int64


        function Parameters(Nx,Ny,beta,spec_nb,rand_seed,is_non_reserve,is_callbacks,is_damier,is_rmax,Rmax=-1)
            new(Nx,Ny,beta,spec_nb,rand_seed,is_non_reserve,is_callbacks,is_damier,is_rmax,Rmax)
        end
    end

# ==============================================================================
# 2 - GRID GRAPH
# ==============================================================================
    struct GridGraph
        Noeuds::Vector{Int}
        Aretes::Vector{Pair{Int,Int}}
        Voisins::Vector{Vector{Int}}
        Arcs::Vector{Pair{Int,Int}}
        NoeudsPeripheriques::Vector{Int}
        NoeudsFictif::Vector{Int}
        AretesFictif::Vector{Pair{Int,Int}}
        VoisinsFictif::Vector{Vector{Int}}
        ArcsFictif::Vector{Pair{Int,Int}}
        alpha::Int64
        NoeudsNoirs::Vector{Int}
        NoeudsBlancs::Vector{Int}

        # Fonction donnant les éléments du graphe de la grille
        function GridGraph(params)

            # Parameters
            Nx        = params.Nx
            Ny        = params.Ny
            beta      = params.beta
            spec_nb   = params.spec_nb
            rand_seed = params.rand_seed

            # Noeuds du graphe associé à la grille de taille Nx*Ny
            Noeuds = 1:Nx*Ny;
            NoeudsFictif = 1:(Nx*Ny+1);

            # Noeud fictif
            alpha = Nx*Ny+1;

            # Liste des voisins pour chaque noeud
            Voisins = Vector{Vector{Int}}()
            for k in Noeuds
                push!(Voisins, Vector{Int}())
            end

            # Arêtes du graphe associé à la grille de taille Nx*Ny
            Aretes = Vector{Pair{Int,Int}}()
            for k in Noeuds
                # arêtes horizontales
                if mod(k,Nx) != 0
                    push!(Aretes,(k=>k+1))
                    push!(Voisins[k], k+1)
                    push!(Voisins[k+1], k)
                end
                # arêtes verticales
                if k <= (Ny-1)*Nx
                    push!(Aretes,(k=>k+Nx))
                    push!(Voisins[k], k+Nx)
                    push!(Voisins[k+Nx], k)
                end
            end

            # Arcs du graphe associé à la grille de taille Nx*Ny
            Arcs = Vector{Pair{Int,Int}}()
            for a in Aretes
            push!(Arcs,a[1]=>a[2])
            push!(Arcs,a[2]=>a[1])
            end

            # Noeuds périphériques du graphe
            NoeudsPeripheriques = Vector{Int64}()
            for k in Noeuds
                if length(Voisins[k]) < 4
                    push!(NoeudsPeripheriques,k)
                end
            end

            # Arêtes du graphe associé à la grille de taille Nx*Ny
            AretesFictif = Vector{Pair{Int,Int}}()
            for k in Noeuds
                # arêtes noeuds fictif avec les noeuds périphériques
                if sum(findall(NoeudsPeripheriques .== k))>0
                    push!(AretesFictif,(alpha=>k))
                end
            end

            # Liste des voisins pour chaque noeud
            VoisinsFictif = Vector{Vector{Int}}()
            for k in NoeudsFictif
                push!(VoisinsFictif, Vector{Int}())
            end

            for k in Noeuds
                # arêtes horizontales
                if mod(k,Nx) != 0
                    push!(VoisinsFictif[k], k+1)
                    push!(VoisinsFictif[k+1], k)
                end
                # arêtes verticales
                if k <= (Ny-1)*Nx
                    push!(VoisinsFictif[k], k+Nx)
                    push!(VoisinsFictif[k+Nx], k)
                end

                # arêtes noeuds fictif avec les noeuds périphériques
                if sum(findall(NoeudsPeripheriques .== k))>0
                    push!(VoisinsFictif[k], alpha)
                    push!(VoisinsFictif[alpha], k)
                end
            end

            # Arcs du graphe associé à la grille de taille Nx*Ny
            ArcsFictif = Vector{Pair{Int,Int}}()
            for a in AretesFictif
                push!(ArcsFictif,a[1]=>a[2])
                push!(ArcsFictif,a[2]=>a[1])
            end
            N_arcsfictif = length(ArcsFictif)

            # Noeuds noirs et blanc du damier
            damier = Damier(params)
            NoeudsNoirs  = damier.Black
            NoeudsBlancs = damier.White

            new(Noeuds,Aretes,Voisins,Arcs,NoeudsPeripheriques,NoeudsFictif,AretesFictif,VoisinsFictif,ArcsFictif,alpha,NoeudsNoirs,NoeudsBlancs)
        end
    end




# ==============================================================================
# 3 - INSTANCE
# ==============================================================================
    struct InstanceFromParams
        N_pu::Int64 #
        N_cf::Int64 #
        N_bd::Int64 #
        PlanningUnits::Array{Int64,1}
        ConservationFeatures::Array{Int64,1}
        Cost::Array{Float64,2}
        Amount::Array{Float64,2}
        Targets::Array{Float64,2}
        Rentability::Array{Float64,2}
        BoundaryLength::Dict{Pair{Int,Int},Int}
        BoundaryCorrection::Array{Float64,2}
        Beta::Float64
        IsLockedOut::Array{Int8,2}


        # Fonction d'habillage de l'instance
        function InstanceFromParams(params,gridgraph)

            # Parameters
            Nx        = params.Nx
            Ny        = params.Ny
            N_pu      = Nx*Ny
            beta      = params.beta
            N_cf      = params.spec_nb
            rand_seed = params.rand_seed

            # Données du graphe de la grille
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
            Cost = zeros(N_noeuds,1)
            IsLockedOut = zeros(N_pu,1)
            for j in Noeuds
                Cost[j] = rand(1:5)
            end

            # La réserve doit recouvrir chaque élément i dans le respect des quantités cibles
            Targets = zeros(N_cf,1)
            for i in ConservationFeatures
                Targets[i] = 0.2*sum(Amount[i,:])
            end

            # Calcul de la rentabilité d'une cellule
            Rentability = zeros(N_noeuds,1)
            for j in Noeuds
                Rentability[j] = sum(Amount[i,j]/Targets[i] for i in ConservationFeatures)/Cost[j]
            end

            # Longueur de chaque arcs
            BoundaryLength = Dict{Pair{Int,Int},Int}()
            for d in Arcs
                BoundaryLength[d] = 1
            end

            # Correction
            BoundaryCorrection = zeros(N_noeuds,1)
            for j in NoeudsPeripheriques
                BoundaryCorrection[j] = 4 - length(Voisins[j])
            end

            new(N_pu,N_cf,N_bd,PlanningUnits,ConservationFeatures,Cost,Amount,Targets,Rentability,BoundaryLength,BoundaryCorrection,beta,IsLockedOut)

        end
    end

    struct InstanceFromFiles
        N_pu::Int64 #
        N_cf::Int64 #
        N_bd::Int64 #
        PlanningUnits::Array{Int64,1}
        ConservationFeatures::Array{Int64,1}
        Cost::Array{Float64,2}
        Amount::Array{Float64,2}
        Targets::Array{Float64,2}
        Rentability::Array{Float64,2}
        BoundaryLength::Dict{Pair{Int,Int},Int}
        BoundaryCorrection::Array{Float64,2}
        Beta::Float64
        IsLockedOut::Array{Int8,2}

        # input data initialisation
        function InstanceFromFiles(pu_fname,cf_fname,bound_fname,beta,targets,gridgraph)

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
            Cost = zeros(N_pu,1)
            IsLockedOut = zeros(N_pu,1)
            #LockedOut = Vector{Int64}()
            for j in 1:N_pu
                Cost[j,1] = pu_data.cost[j]
                IsLockedOut[j,1] = pu_data.is_locked_out[j]
            end

            # Amount and targets of conservation features
            Amount = zeros(N_cf,N_pu)
            Targets = zeros(N_cf,1)
            for i in 1:N_cf
                Amount[i,1:N_pu] = cf_data[:,i+1]
                Targets[i,1] = targets[i] * sum(Amount[i,1:N_pu] .* (1 .- IsLockedOut))
            end

            # Calcul de la rentabilité d'une cellule
            Rentability = zeros(N_noeuds,1)
            for j in Noeuds
                Rentability[j] = sum(Amount[i,j]/Targets[i] for i in ConservationFeatures)/Cost[j]
            end

            # Boundary length of vertices between two nodes
            BoundaryLength = Dict{Pair{Int,Int},Int}()
            for k in 1:N_bd
                BoundaryLength[bound_data.id1[k]=>bound_data.id2[k]] = bound_data.boundary[k]
                BoundaryLength[bound_data.id2[k]=>bound_data.id1[k]] = bound_data.boundary[k]
            end

            # Correction
            BoundaryCorrection = zeros(N_pu,1)
            for j in NoeudsPeripheriques
                BoundaryCorrection[j] = 4 - length(Voisins[j])
            end

            # constructor
            new(N_pu,N_cf,N_bd,PlanningUnits,ConservationFeatures,Cost,Amount,Targets,Rentability,BoundaryLength,BoundaryCorrection,beta,IsLockedOut)

        end
    end

# ==============================================================================
# 4 - DAMIER
# ==============================================================================
    struct Damier
        IsBlack::Array{Int64,1}
        Black::Vector{Int64}
        White::Vector{Int64}

        # Fonction générant la liste des noeuds noirs et blancs du damier de la grille
        function Damier(params)

            # Parameters
            Nx = params.Nx
            Ny = params.Ny

            # Initialisation
            IsBlack = Array{Int64,1}(undef,Nx*Ny)
            Black   = Vector{Int64}()
            White   = Vector{Int64}()

            # Damier
            num_ligne = 0
            for k in 1:(Nx*Ny)
                if mod(k,Nx) == 0
                    num_colonne = Nx
                else
                    num_colonne = mod(k,Nx)
                end

                if mod(k,Nx) == 1
                    num_ligne = num_ligne + 1
                end

                IsBlack[k] = 0
                if mod(num_ligne,2) == mod(num_colonne,2)
                    IsBlack[k] = 1
                    push!(Black,k)
                else
                    push!(White,k)
                end
            end

            new(IsBlack, Black, White)
        end
    end
