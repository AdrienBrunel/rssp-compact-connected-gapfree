# ==============================================================================
# 1 - PRINT INFO
# ==============================================================================


    # Fonction générant la liste des noeuds noirs et blancs du damier de la grille
    function print_grid_graph(gridgraph)

        # Données du graphe de la grille
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

        # Affichage
        println("Listes de $(length(Noeuds)) noeuds : $Noeuds")
        println("Liste des $(length(Aretes)) arêtes: $Aretes")
        println("Listes des voisins : $Voisins")
        println("Liste des $(length(Arcs)) arcs: $Arcs")
        println("Liste des $(length(NoeudsPeripheriques)) noeuds périphériques: $NoeudsPeripheriques")
        println("Listes de $(length(NoeudsFictif)) noeuds : $NoeudsFictif")
        println("Liste des $(length(AretesFictif)) arêtes: $AretesFictif")
        println("Listes des voisins fictifs : $VoisinsFictif")
        println("Liste des $(length(ArcsFictif)) arcs: $ArcsFictif")
        println("Le noeud fictif est le noeud $alpha")

        return
    end


    # Fonction générant la liste des noeuds noirs et blancs du damier de la grille
    function print_instance(instance,gridgraph)

        # Données de l'instance
        ConservationFeatures = instance.ConservationFeatures
        Amount               = instance.Amount
        Cost                 = instance.Cost
        Targets              = instance.Targets
        Rentability          = instance.Rentability
        BoundaryLength       = instance.BoundaryLength
        BoundaryCorrection   = instance.BoundaryCorrection

        # Données du graphe de la grille
        Noeuds = gridgraph.Noeuds
        Arcs   = gridgraph.Arcs

        # Affichage
        println(Amount[1,:])
        println(Amount[2,:])
        println(Amount[3,:])
        println(Cost)
        println(Targets)
        println(Rentability)
        for j in Noeuds
            if BoundaryCorrection[j] > 0
                println("Le noeud $j a une correction de sa Boundary Length de $(BoundaryCorrection[j])")
            end
        end
        println(sum(BoundaryLength[d] for d in Arcs))

        return
    end


    # Fonction affichant la liste des noeuds noirs et blancs du damier de la grille
    function print_damier(damier)

        # Récupération des infos du damier
        IsBlack = damier.IsBlack
        Black   = damier.Black
        White   = damier.White

        # Affichage
        println("Listes des noeuds noirs : $Black")
        println("Listes des noeuds blancs : $White")

        return
    end


# ==============================================================================
# 2 - PLAN COUPANTS
# ==============================================================================

    # Fonction donnant la liste des composantes connexes
    function connected_components(neighbors::Vector{Vector{Int64}}, is_reserve::Vector{Int64})
        # get the support of current solution
        R = findall(is_reserve .>= 1.0e-4)
        nv = length(is_reserve)
        components = Vector{Vector{Int64}}()
        marked = falses(nv)
        for u in R

            if marked[u] continue end
            component = Vector{Int64}()
            queue = Vector{Int64}()
            push!(queue, u)
            marked[u] = true

            # breadth first search starting from u provides the component including u
            while !isempty(queue)
                v = popfirst!(queue)
                push!(component, v)
                for w in neighbors[v]
                    if is_reserve[w] >= 1.0e-4 && !marked[w]
                        push!(queue, w)
                        marked[w] = true
                    end
                end
            end
            push!(components, component)
        end
        return components
    end

    # Compute the distance between each pair of vertices
    function shortest_distances(nv::Int64, neighbors::Vector{Vector{Int64}})
        # initialize pairwise distances
        dmin = nv * ones(Int, nv, nv)

        for u in 1:nv
            # compute distance from u to every other vertex
            dmin[u,u] = 0
            marked = falses(nv)
            marked[u] = true
            queue = Vector{Int64}()
            push!(queue, u)
            # every arc has length 1 so a breadth-first search is enough
            while !isempty(queue)
                v = popfirst!(queue)
                for w in neighbors[v]
                    if !marked[w]
                        dmin[u,w] = dmin[u,v] + 1
                        push!(queue, w)
                        marked[w] = true
                    end
                end
            end
        end
        return dmin
    end

    # Fonction calculant le plan coupant à introduire
    function get_cut(nv::Int64, idx_C1::Int64, components::Vector{Vector{Int64}}, dmin::Matrix{Int64}, neighbors::Vector{Vector{Int64}})

        C1 = components[idx_C1]

        # compute belonging indicators for sets C1 and C2
        is_in_C1 = falses(nv)
        for u in C1
            is_in_C1[u] = true
        end

        # compute the frontier of sets C1 and C2
        δC1 = Vector{Int64}()
        for u in C1
            for v in neighbors[u]
                if !is_in_C1[v]
                    push!(δC1,u)
                    #push!(v, δC1)
                    break
                end
            end
        end

        # compute shortest distance from C1 and C2 to every other vertex
        dC1 = minimum(dmin[δC1, :], dims=1)
        @assert(length(dC1) == nv)
        for u in C1
            dC1[u] = 0
        end

        dC1_min = nv
        idx_C2 = 0
        for i in 1:length(components)
            if i==idx_C1
                continue
            end
            Ci = components[i]
            dC1Ci = minimum(dC1[Ci])
            if dC1Ci < dC1_min
                dC1_min = dC1Ci
                idx_C2 = i
            end
        end
        C2 = components[idx_C2]
        println("C2 = $(C2)")

        dC1C2 = minimum(dC1[C2])
        @assert(dC1C2 > 0)
        println("dC1C2 = $(dC1C2)")

        # define the cuts
        S1 = C1
        S2 = Vector{Int64}()
        Sbar = Vector{Int64}()
        for u in 1:nv
            if (dC1[u] == 0)
                continue
            elseif (dC1[u] < dC1C2)
                push!(Sbar,u)
            else
                push!(S2,u)
            end
        end
        return S1, S2, Sbar, dC1C2
    end


# ==============================================================================
# 3 - LAYOUT GRAPHS
# ==============================================================================

    # Fonction permettant la visualisation propre de la grille sous forme de graphe
    function grid_layout(g)
        locs_x = Array{Int64,1}(undef,Nx*Ny)
        locs_y = Array{Int64,1}(undef,Nx*Ny)

        num_ligne = 0
        for k in 1:(Nx*Ny)
            if mod(k,Nx) == 0
                num_colonne = Nx
            else
                num_colonne = mod(k,Nx)
            end
            locs_x[k] = num_colonne

            if mod(k,Nx) == 1
                num_ligne = num_ligne + 1
            end
            locs_y[k] = num_ligne
        end
        return locs_x, locs_y
    end


    # Fonction permettant la visualisation propre de la grille sous forme de graphe avec le noeud fictif
    function grid_layout_fictif(g)
        locs_x = Array{Int64,1}(undef,Nx*Ny+1)
        locs_y = Array{Int64,1}(undef,Nx*Ny+1)

        num_ligne = 0
        for k in 1:(Nx*Ny)
            if mod(k,Nx) == 0
                num_colonne = Nx
            else
                num_colonne = mod(k,Nx)
            end
            locs_x[k] = num_colonne

            if mod(k,Nx) == 1
                num_ligne = num_ligne + 1
            end
            locs_y[k] = num_ligne
        end
        locs_x[Nx*Ny+1]=0
        locs_y[Nx*Ny+1]=0
        return locs_x, locs_y
    end



# ==============================================================================
# 5 - LECTURE DES RESULTATS
# ==============================================================================

    function read_info_reserve(model,gridgraph)

        # Lecture données du graphe de la grille
        Arcs   = gridgraph.Arcs
        Noeuds = gridgraph.Noeuds

        # variable de sélection des noeuds de la reserve
        x_opt = round.(Int,value.(model[:x]).data)
        println("La réserve est composée des noeuds $(Noeuds[x_opt[Noeuds] .== 1])")

        # variable de linéarisation
        z_tmp = round.(Int,value.(model[:z]).data)
        z_opt = Dict{Pair{Int,Int},Int}()
        for d in Arcs
            z_opt[d] = 0
            if sum(findall(d .== Arcs[z_tmp .==1]))>0
                z_opt[d] = 1
            end
        end

        return x_opt, z_opt
    end


    # Lecture des infos
    function read_info_graph_reserve(model,gridgraph,is_beta)

        # Lecture données du graphe de la grille
        Arcs   = gridgraph.Arcs
        Noeuds = gridgraph.Noeuds

        # variable de sélection des noeuds de la reserve
        x_opt = round.(Int,value.(model[:x]).data)
        println("La réserve est composée des noeuds $(Noeuds[x_opt[Noeuds] .== 1])")

        # variable de sélection de la racine de l'arbre
        r_opt = round.(Int,value.(model[:r]).data)
        println("Le noeud $(Noeuds[r_opt .== 1][1]) est la racine de l'arbre couvrant de la réserve")

        # variable de linéarisation
        if is_beta
            z_tmp = round.(Int,value.(model[:z]).data)
            z_opt = Dict{Pair{Int,Int},Int}()
            for d in Arcs
                z_opt[d] = 0
                if sum(findall(d .== Arcs[z_tmp .==1]))>0
                    z_opt[d] = 1
                end
            end
        else
            z_opt = -1
        end

        # variable de sélection d'arcs
        u_tmp = round.(Int,value.(model[:u]).data)
        u_opt = Dict{Pair{Int,Int},Int}()
        for d in Arcs
            u_opt[d] = 0
            if sum(findall(d .== Arcs[u_tmp .==1]))>0
                u_opt[d] = 1
            end
        end
        println("Les arcs constituant l'arbre couvrant de la réserve sont $(Arcs[u_tmp .== 1])")

        return x_opt,z_opt,u_opt,r_opt
    end


    # Lecture des infos
    function read_info_graph_nonreserve(model,gridgraph)

        # Lecture données du graphe de la grille
        Arcs         = gridgraph.Arcs
        Noeuds       = gridgraph.Noeuds
        NoeudsFictif = gridgraph.NoeudsFictif
        ArcsFictif   = gridgraph.ArcsFictif

        # variable de sélection des noeuds de la reserve
        x_opt = round.(Int,value.(model[:x]).data)
        push!(x_opt, 0)
        println("La non-réserve est composée des noeuds $(NoeudsFictif[x_opt[NoeudsFictif] .== 0])")

        # variable de sélection de la racine de l'arbre
        println("Le noeud fictif $(alpha) est la racine de l'arbre couvrant de la non-réserve")

        # variable de sélection d'arcs
        v_tmp = round.(Int,value.(model[:v]).data)
        v_opt = Dict{Pair{Int,Int},Int}()
        for d in union(Arcs,ArcsFictif)
            v_opt[d] = 0
            if sum(findall(d .== union(Arcs,ArcsFictif)[v_tmp .==1]))>0
                v_opt[d] = 1
            end
        end
        println("Les arcs constituant l'arbre couvrant de la non-réserve sont $(union(Arcs,ArcsFictif)[v_tmp .== 1])")

        return x_opt,v_opt
    end


    # Affichage des infos
    function print_info_reserve(x_opt,instance,gridgraph)

        # Lecture données du graphe de la grille
        Noeuds              = gridgraph.Noeuds
        Voisins             = gridgraph.Voisins
        NoeudsPeripheriques = gridgraph.NoeudsPeripheriques
        N_noeuds = length(Noeuds)

        # Lecture données de l'instance
        ConservationFeatures = instance.ConservationFeatures
        Amount               = instance.Amount
        Cost                 = instance.Cost
        Targets              = instance.Targets
        BoundaryLength       = instance.BoundaryLength
        BoundaryCorrection   = instance.BoundaryCorrection

        # variable de sélection des noeuds de la reserve
        println("La réserve est composée des noeuds $(Noeuds[x_opt[Noeuds] .== 1])")
        println("La réserve coûte $(sum(Cost[x_opt[Noeuds] .== 1]))")
        for i in ConservationFeatures
            println("La réserve contient $(sum(Amount[i,x_opt[Noeuds] .== 1])) de l'espèce $i dont la cible est $(Targets[i])")
        end

        # Elements de la reserve
        Perimetre=0
        for j in Noeuds
            if x_opt[j] == 1
                Perimetre = Perimetre + sum(BoundaryLength[j=>k]*(1-x_opt[k]) for k in Voisins[j])+BoundaryCorrection[j]*x_opt[j]
            end
        end
        Cout = sum(Cost[x_opt[Noeuds] .== 1])
        Score = Cout + beta*Perimetre
        println("Le périmètre de la réserve est $Perimetre")
        println("Le score de la réserve est : Score = Cout + beta x Perimetre = $(Cout) + $(beta) x $(Perimetre) = $(Score)")

        # lecture des données du graphe de la grille
        dmin = shortest_distances(N_noeuds,Voisins)
        Reserve = Noeuds[findall(x_opt .== 1)]
        Rayon = 0
        for j in Reserve
            tmp = maximum(dmin[j,Reserve]/2)
            if tmp > Rayon
                Rayon = tmp
            end
        end
        println("- outside radius of the reserve: $(Rayon)")

        # compute the inside radius of the reserve
        # - for each node of the reserve, compute its neighbors in the reserve
        Voisins_InReserve = Vector{Vector{Int}}()
        for i in 1:length(Reserve)
            push!(Voisins_InReserve, Vector{Int}())
            n1 = Reserve[i]
            for n2 in  Reserve
                if n2 ∈ Voisins[n1]
                    push!(Voisins_InReserve[i], findfirst(Reserve .== n2))
                end
            end
        end
        # - get minimum pairwise distances
        dmin_in_reserve = shortest_distances(length(Reserve),Voisins_InReserve)
        # - the center of reserve is the node whose maximum distance to other nodes is smallest
        max_dmin = maximum(dmin_in_reserve, dims = 2)
        Rayon_InReserve, center = findmin(max_dmin)
        println("- inside radius of the reserve: $(Rayon_InReserve)")

        # verify that the reserve is connected
        Reserve_Components = connected_components(Voisins,x_opt)
        if length(Reserve_Components) == 1
            println("- the reserve is connected")
        else
            println("- the reserve has $(length(Reserve_Components)) connected components")
        end

        # detect the number of holes in the reserve
        VoisinsFictif  = gridgraph.VoisinsFictif
        # - add an artificial node in the solution if not already there
        if length(x_opt) == N_noeuds
            push!(x_opt, 0)
        end
        # - get the number of connected components in the complementary of the reserve
        NonReserve_Components = connected_components(VoisinsFictif,1 .- x_opt)
        if length(NonReserve_Components) == 1
            println("- the reserve has no hole")
        else
            println("- the reserve has $(length(NonReserve_Components)-1) holes")
        end

        return Perimetre,Cout,Score,Rayon,Rayon_InReserve
    end


# ==============================================================================
# 5 - VISUALISATION
# ==============================================================================


    function visualisation_reserve(x_opt,dir,gridgraph)

        # lecture des données du graphe de la grille
        Noeuds   = gridgraph.Noeuds
        Arcs     = gridgraph.Arcs
        N_noeuds = length(Noeuds)

        # création de l'objet graph de julia
        G = DiGraph(N_noeuds)
        nodelabels = collect(1:N_noeuds)

        for a in Arcs
            add_edge!(G,a[1],a[2])
        end

        nodefillc = []
        for i in Noeuds
            if x_opt[i]==1
                push!(nodefillc, RGBA(0.0,1.0,0.0,0.5))
            else
                push!(nodefillc, RGBA(0.0,0.0,1.0,0.5))
            end
        end

        draw(PNG(dir, 36cm, 21cm),
        gplot(G,nodelabel=nodelabels,nodefillc=nodefillc,layout=grid_layout,arrowlengthfrac=0.05,edgelinewidth=1))
        return

    end



    function visualisation_reserve_graph(x_opt,u_opt,r_opt,dir,gridgraph)

        # lecture des données du graphe de la grille
        Noeuds   = gridgraph.Noeuds
        Arcs     = gridgraph.Arcs
        N_noeuds = length(Noeuds)

        # création de l'objet graph de julia
        G = DiGraph(N_noeuds)
        nodelabels = collect(1:N_noeuds)

        for a in Arcs
            add_edge!(G,a[1],a[2])
        end

        nodefillc = []
        for i in Noeuds
            if x_opt[i]==1
                if r_opt[i]==1
                    push!(nodefillc, RGBA(0.0,1.0,0.0,1.0))
                else
                    push!(nodefillc, RGBA(0.0,1.0,0.0,0.5))
                end
            else
                push!(nodefillc, RGBA(0.0,0.0,1.0,0.5))
            end
        end

        edgestrokec = []
        for a in edges(G)
            if u_opt[(src(a)=>dst(a))] == 1
                push!(edgestrokec,RGBA(1.0,0.0,0.0,1.0))
            else
                push!(edgestrokec,RGBA(0.0,0.0,0.0,0.0))
            end
        end

        draw(PNG(dir, 36cm, 21cm),
        gplot(G,nodelabel=nodelabels,nodefillc=nodefillc,layout=grid_layout,edgestrokec=edgestrokec,arrowlengthfrac=0.05,edgelinewidth=1))
        return

    end




    function visualisation_reserve_nonreserve_graph(x_opt,u_opt,r_opt,v_opt,dir,gridgraph)

        # lecture des données du graphe de la grille
        NoeudsFictif   = gridgraph.NoeudsFictif
        Arcs           = gridgraph.Arcs
        ArcsFictif     = gridgraph.ArcsFictif
        N_noeudsfictif = length(NoeudsFictif)

        # création de l'objet graph de julia
        G = DiGraph(N_noeudsfictif)
        nodelabels = collect(1:N_noeudsfictif)

        for a in union(Arcs,ArcsFictif)
            add_edge!(G,a[1],a[2])
        end

        nodefillc = []
        for i in NoeudsFictif
            if x_opt[i]==1
                if r_opt[i]==1
                    push!(nodefillc, RGBA(0.0,1.0,0.0,1.0))
                else
                    push!(nodefillc, RGBA(0.0,1.0,0.0,0.5))
                end
            else
                push!(nodefillc, RGBA(0.0,0.0,1.0,0.5))
            end
        end

        edgestrokec = []
        for a in edges(G)
            if sum(findall(Arcs .== (src(a)=>dst(a))))>0
                if u_opt[(src(a)=>dst(a))] == 1
                    push!(edgestrokec,RGBA(1.0,0.0,0.0,1.0))
                elseif v_opt[(src(a)=>dst(a))] == 1
                    push!(edgestrokec,RGBA(1.0,0.5,0.0,1.0))
                else
                    push!(edgestrokec,RGBA(0.0,0.0,1.0,0.0))
                end
            else
                if v_opt[(src(a)=>dst(a))] == 1
                    push!(edgestrokec,RGBA(1.0,0.5,0.0,1.0))
                else
                    push!(edgestrokec,RGBA(0.0,0.0,1.0,0.0))
                end
            end
        end

        draw(PNG(dir, 36cm, 21cm),
        gplot(G,nodelabel=nodelabels,nodefillc=nodefillc,layout=grid_layout_fictif,edgestrokec=edgestrokec,arrowlengthfrac=0.05,edgelinewidth=1))
        return

    end
