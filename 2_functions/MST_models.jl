# ==============================================================================
# 0 - DEV
# ==============================================================================
    function ReserveSiteSelection_SpatialConstraints(instance, gridgraph, params, is_beta, is_non_reserve, is_callbacks, is_damier, is_rmax)

        ## LECTURE DES DONNEES -------------------------------------------------
        Noeuds               = gridgraph.Noeuds
        Voisins              = gridgraph.Voisins
        Arcs                 = gridgraph.Arcs
        NoeudsPeripheriques  = gridgraph.NoeudsPeripheriques
        N_noeuds             = length(Noeuds)
        ConservationFeatures = instance.ConservationFeatures
        Amount               = instance.Amount
        Cost                 = instance.Cost
        Targets              = instance.Targets
        BoundaryLength       = instance.BoundaryLength
        BoundaryCorrection   = instance.BoundaryCorrection
        Beta                 = instance.Beta
        Rmax = params.Rmax

        if is_damier
            NoeudsNoirs  = gridgraph.NoeudsNoirs
        end

        if is_non_reserve
            VoisinsFictif  = gridgraph.VoisinsFictif
            ArcsFictif     = gridgraph.ArcsFictif
            alpha          = gridgraph.alpha
        end

        ## VARIABLE DE DECISION ------------------------------------------------
        # x : variable de sélection du noeuds j dans le graphe de la réserve
        # z : variable de linéarisation du périmètre quadratique
        # u : variable de sélection de l'arc dans le graphe de la réserve
        # r : variable de sélection du noeud racine de l'arbre couvrant
        # v : variable de sélection de l'arc dans le graphe de la non-réserve
        # f : variable de sélection du flux dans le graphe de la réserve
        # g : variable de sélection du flux dans le graphe de la non-réserve

        m = Model(Gurobi.Optimizer)
		set_optimizer_attribute(m, "TimeLimit", 1000)
		set_optimizer_attribute(m, "LogFile", "$(res_dir)/gurobi_log.txt")

        NoeudsInterieurs = setdiff(Noeuds,NoeudsPeripheriques)
        if is_damier
            Commodites_nr = intersect(NoeudsNoirs, NoeudsInterieurs)
        else
            Commodites_nr = NoeudsInterieurs
        end
        Commodites = Noeuds
		if is_beta
        	@variable(m, z[Arcs], Bin);
		end
        @variable(m, u[Arcs], Bin);
        @variable(m, r[Noeuds], Bin);
        if is_non_reserve
            @variable(m, v[union(Arcs,ArcsFictif)], Bin);
        end
        @variable(m, x[Noeuds], Bin);

        # define the flow variables for connexity of reserve and non-reserve
        # - get the nodes that will be concerned by each reserve commodity
        dmin = shortest_distances(N_noeuds,Voisins)
        NoeudsProche  = Vector{Vector{Int}}()
        for k in Noeuds
            push!(NoeudsProche, Vector{Int}())
            for j in Noeuds
                if dmin[k,j]<=Rmax
                    push!(NoeudsProche[k],j)
                end
            end
        end
        @variable(m, f[k in Commodites, i in NoeudsProche[k], j in Voisins[i];dmin[k,i]<=Rmax-1] >= 0);
        # - take care of non-reserve flow variables
        if is_non_reserve
            NoeudsProche_nr  = Vector{Vector{Int}}()
            NoeudsFrontiere_nr = Vector{Vector{Int}}()
            for k in Noeuds
                push!(NoeudsProche_nr, Vector{Int}())
                push!(NoeudsFrontiere_nr, Vector{Int}())
                for j in Noeuds
                    # Rmax is used below but 2Rmax-2 should be used for theoretical proof
                    if dmin[k,j]<= Rmax
                        push!(NoeudsProche_nr[k],j)
                        if dmin[k,j] == Rmax
                            push!(NoeudsFrontiere_nr[k],j)
                        end
                    end
                end
            end
            @variable(m, g[k in Commodites_nr, i in NoeudsProche_nr[k], j in Voisins[i]] >= 0);
        end

		## OBJECTIF ------------------------------------------------------------
		if is_beta
	    	@objective(m, Min, sum(Cost[j]*x[j] for j in Noeuds) + beta*sum(BoundaryLength[d]*(x[d[1]]-z[d]) for d in Arcs) + Beta*sum(BoundaryCorrection[j]*x[j] for j in NoeudsPeripheriques))
		else
			@objective(m, Min, sum(Cost[j]*x[j] for j in Noeuds))
		end

        ## CONTRAINTES ---------------------------------------------------------

		# Reserve site selection
        @constraint(m, cibles[i in ConservationFeatures], sum(Amount[i,j]*x[j] for j in Noeuds) >= Targets[i])

        # Linearization of the perimeter term in the objective
		if is_beta
	        @constraint(m, linzi[d in Arcs],  z[d] - x[d[1]] <= 0)
        	@constraint(m, linzk[d in Arcs],  z[d] - x[d[2]] <= 0)
		end

        # The nodes indicated by x and the edges indicated by u must describe a reverse directed tree rooted at the node indicated by r
        # @constraint(m, arc_unique[a in Aretes],  u[a[1]=>a[2]] + u[a[2]=>a[1]] <= 1) # useless/redundant
        @constraint(m, nb_arcs, sum(u[d] for d in Arcs) == sum(x[j] for j in Noeuds)-1)
        @constraint(m, arc_reserve_1[d in Arcs],  u[d] <= x[d[1]])
        @constraint(m, arc_reserve_2[d in Arcs],  u[d] <= x[d[2]])
        @constraint(m, racine_reserve[j in Noeuds],  r[j] <= x[j])
        @constraint(m, racine_unique, sum(r[j] for j in Noeuds) <= 1)
        # In a reverse directed tree, each node that is not the root has exactly one out-arc
        @constraint(m, visu_r[j in Noeuds], sum(u[j=>i] for i in Voisins[j]) ==  x[j]-r[j])

        # Arbre couvrant de la non-reserve
        if is_non_reserve
            # @constraint(m, arc_unique_nr[a in union(Aretes,AretesFictif)],  v[a[1]=>a[2]] + v[a[2]=>a[1]] <= 1) # useless/redundant
            @constraint(m, nb_arcs_nr, sum(v[d] for d in union(Arcs,ArcsFictif)) == sum(1-x[j] for j in Noeuds))
            @constraint(m, arc_reserve_1_nr[i in NoeudsInterieurs,j in Voisins[i]],  v[i=>j] <= 1-x[i])
            @constraint(m, arc_reserve_2_nr[i in NoeudsInterieurs,j in Voisins[i]],  v[i=>j] <= 1-x[j])
            @constraint(m, one_out_arc_int[j in NoeudsInterieurs], sum(v[j=>i] for i in Voisins[j]) ==  1-x[j])
            @constraint(m, no_out_arc_border[j in NoeudsPeripheriques], sum(v[j=>i] for i in Voisins[j]) == 0)
            @constraint(m, no_arc_from_alpha[d in ArcsFictif; d[1] == alpha], v[d] == 0)
            @constraint(m, all_arcs_to_alpha[d in ArcsFictif; d[2] == alpha], v[d] == 1 - x[d[1]])
        end

        # avoid isolated nodes in the reserve and in the non-reserve (former Réduction damier but for all nodes)
        @constraint(m, no_isolated[j in Noeuds], x[j] <= sum(x[i] for i in Voisins[j]))
        @constraint(m, no_isolated_nr[j in NoeudsInterieurs], 1-x[j] <= sum(1-x[i] for i in Voisins[j]))

        if is_rmax
            @constraint(m, source_rmax[j in Noeuds, i in Noeuds;dmin[i,j]>Rmax], x[i] <=  1-r[j])
        end

        if is_damier
            Commodites_restante = copy(NoeudsNoirs)
        else
            Commodites_restante = copy(Commodites)
        end
        Commodites_restante_nr = copy(Commodites_nr)
        Commodites_choisies = []
        Commodites_choisies_nr = []

        function my_callback_function(cb_data)
            # stop the callback if the solution is non-integer
            x_val = zeros(Int,N_noeuds)
            # x_val = Dict{Int,Int}()
            for j in Noeuds
                val = callback_value(cb_data, x[j])
                if abs(val-round(Int,val)) < 1e-6
                    x_val[j] = round(Int,val)
                else
                    return
                end
                val = callback_value(cb_data, r[j])
                if abs(val-round(Int,val)) >= 1e-6
                    return
                end
            end
            for d in Arcs
                val = callback_value(cb_data, u[d])
                if abs(val-round(Int,val)) >= 1e-6
                    return
                end
                if is_non_reserve
                    val = callback_value(cb_data, v[d])
                    if abs(val-round(Int,val)) >= 1e-6
                        return
                    end
                end
            end
            if is_non_reserve
                push!(x_val, 0)
            end

            # get the connected components of the reserve
            Reserve_Components = connected_components(Voisins,x_val)

            # select one node per component to add it as a commodity
            Commodities_to_add = []
            if length(Reserve_Components) > 1
                println("cb: $(length(Reserve_Components)) connected components in the reserve : $(Reserve_Components)")

                for c in 1:length(Reserve_Components)

                    # choix de la commodité
                    Component = intersect(Reserve_Components[c],Commodites_restante)
                    # make sure that a commodity is not taken twice
                    if isempty(Component)
                        println("cb: every commodity of the component has already been added")
                        continue
                    end 
                    k = Component[findmax(instance.Rentability[Component])[2][1]]
                    push!(Commodities_to_add, k)
                    push!(Commodites_choisies, k)
                    setdiff!(Commodites_restante, k)
                    println("cb: chosen commodity $(k)")
                end
            end

            # if the reserve is connected but needs to have a maximum radius compute pairwise distances in the reserve
            if is_rmax && (length(Reserve_Components) == 1)
                Reserve = Noeuds[findall(x_val .== 1)]
                
                # for each node of the reserve, compute its neighbors in the reserve
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

                # get minimum pairwise distances
                dmin_in_reserve = shortest_distances(length(Reserve),Voisins_InReserve)
                # get true center of reserve: it is the node whose maximum distance from other nodes is smallest
                max_dmin = maximum(dmin_in_reserve, dims = 2)
                (radius, center) = findmin(max_dmin)
                center = center[1]
                # add the commodities of nodes that are too far away from true center
                # WARNING: possible bug here if all the commodities far from the center are already added
                if radius > Rmax
                    println("cb: center is $(Reserve[center]) and radius is $radius")
                    for k in Reserve[findall(dmin_in_reserve[center,:] .== radius)]
                        if k ∈ Commodites_choisies
                            continue
                        else
                            push!(Commodites_choisies, k)
                            setdiff!(Commodites_restante, k)
                            println("cb: node $k too far away, add commodity")
                            push!(Commodities_to_add, k)
                            break
                        end
                    end
                end 
            end

            # add the lazy cuts for connectivity and size oif the reserve
            for k in Commodities_to_add
                # Ajout des contraintes de flux de la réserve
                for i in NoeudsProche[k]
                    if dmin[k,i]==Rmax
                        # if i is exactly Rmax away from k, no flow can go out from i, so the in-flow can be equal to one iff it is the root
                        conservation_flux_1 = @build_constraint(sum(f[k,j,i] for j in Voisins[i] if dmin[k,j] <= Rmax-1) <= r[i])
                        MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1)
                        continue
                    end
                    for j in Voisins[i]
                        # flow can only go through selected arcs
                        flux_arc_1 = @build_constraint(f[k,i,j] <= u[i=>j])
                        MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_1)
                        # every flow is zero if the commodity is not in the reserve
                        flux_arc_2 = @build_constraint(f[k,i,j] <= x[k])
                        MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_2)
                    end
                    if i==k
                        continue
                    end
                    # flow is conserved at every node but the root, where in-flow can be equal to one (if the flow from k is also equal to one)
                    conservation_flux_1 = @build_constraint(sum(f[k,j,i] for j in Voisins[i] if dmin[k,j] <= Rmax-1) - sum(f[k,i,j] for j in Voisins[i]) <= r[i])
                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1)
                    conservation_flux_2 = @build_constraint(sum(f[k,j,i] for j in Voisins[i] if dmin[k,j] <= Rmax-1) - sum(f[k,i,j] for j in Voisins[i]) >= 0)
                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_2)
                end
                # the in-flow in k is zero
                source_k_in = @build_constraint(sum(f[k,j,k] for j in Voisins[k]) == 0)
                MOI.submit(m, MOI.LazyConstraint(cb_data), source_k_in)
                # 1 value of flow goes out of k iff k is in the reserve and it is not the root of the tree
                source_k_out = @build_constraint(sum(f[k,k,j] for j in Voisins[k]) == x[k]-r[k])
                MOI.submit(m, MOI.LazyConstraint(cb_data), source_k_out)

                if is_rmax
                    # no more than Rmax positive flows for the commodity
                    rmax_inreserve = @build_constraint(sum(f[k, i, j] for i in NoeudsProche[k] for j in Voisins[i] if dmin[k,i] <= Rmax-1) <= Rmax)
                    MOI.submit(m, MOI.LazyConstraint(cb_data), rmax_inreserve)
                end
            end

            # if no lazy cut has been added for the reserve, look at the connectivity of nonreserve
            if isempty(Commodities_to_add) && is_non_reserve
                # Lazy cuts pour la connectivité de la non-réserve
                NonReserve_Components = connected_components(VoisinsFictif,1 .- x_val)
                if length(NonReserve_Components) > 1
                    println("cb: $(length(NonReserve_Components)) connected components in non-reserve: $(NonReserve_Components)")

                    for c in 1:length(NonReserve_Components)
                        if !(sum(findall(NonReserve_Components[c] .== alpha))>1)
                            Component = intersect(NonReserve_Components[c],Commodites_restante_nr)
                            # make sure that a commodity is not taken twice
                            if isempty(Component)
                                println("cb: every commodity of the component has already been added")
                                continue
                            end 
                            # choose commodity
                            k = Component[findmin(instance.Rentability[Component])[2][1]]
                            push!(Commodites_choisies_nr, k)
                            setdiff!(Commodites_restante_nr, k)
                            println("cb: chosen commodity $(k)")

                            # Ajout des contraintes de flux de la non-réserve
                            for i in NoeudsProche_nr[k]
                                for j in Voisins[i]
                                    flux_arc_1_nr = @build_constraint(g[k,i,j] <= v[i=>j])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_1_nr)
                                    flux_arc_1_nr = @build_constraint(g[k,i,j] <= v[i=>j])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_1_nr)
                                    flux_arc_2_nr = @build_constraint(g[k,i,j] <= 1-x[k])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_2_nr)
                                end
                            end
                            # no flow goes in k
                            source_k_nr_in = @build_constraint(sum(g[k,j,k] for j in Voisins[k]) == 0)
                            MOI.submit(m, MOI.LazyConstraint(cb_data), source_k_nr_in)
                            # 1 flow goes out of k if it is not in the reserve
                            source_k_nr_out = @build_constraint(sum(g[k,k,j] for j in Voisins[k]) == 1 - x[k])
                            MOI.submit(m, MOI.LazyConstraint(cb_data), source_k_nr_out)

                            # no flow goes out of alpha
                            Frontiere = intersect(NoeudsPeripheriques, NoeudsProche_nr[k])
                            union!(Frontiere, NoeudsFrontiere_nr[k])
                            Interieur = setdiff(NoeudsProche_nr[k], Frontiere)
                            sink_k_nr_out = @build_constraint(sum(g[k,i,j] for i in Frontiere for j in Voisins[i]) == 0)
                            MOI.submit(m, MOI.LazyConstraint(cb_data),sink_k_nr_out)
                            # 1 flow goes to the frontier if k is not in the reserve
                            sink_k_nr_in = @build_constraint(sum(g[k,j,i] for i in Frontiere for j in Voisins[i] if j in Interieur) == 1 - x[k])
                            MOI.submit(m, MOI.LazyConstraint(cb_data), sink_k_nr_in)

                            # at all other nodes flow must be conserved
                            for i in Interieur
                                if i != k
                                    conservation_flux_nr = @build_constraint(sum(g[k, i, j] for j in Voisins[i]) - sum(g[k,j, i] for j in Voisins[i]) == 0)
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_nr)
                                end
                            end
                            #@constraint(m, conservation_flux_2_nr[k in Noeuds,  n in Noeuds;n!=k], sum(g[k,j=>n] for j in VoisinsFictif[n]) - sum(g[k,n=>j] for j in VoisinsFictif[n]) == 0)
                        end
                    end
                end
            end
        end
        # Contrainte sur les variables de flux
        if is_callbacks
            MOI.set(m, MOI.LazyConstraintCallback(), my_callback_function)
        else
            # # WARNING: this is not compatible with Rmax
            @constraint(m, flux_arc_1[k in Commodites, d in Arcs], f[k,d] <= u[d])
            @constraint(m, flux_arc_2[k in Commodites, d in Arcs], f[k,d] <= x[k])
            @constraint(m, conservation_flux_1[k in Commodites, n in Noeuds;n!=k], sum(f[k,j=>n] for j in Voisins[n]) - sum(f[k,n=>j] for j in Voisins[n]) <= r[n])
            @constraint(m, conservation_flux_2[k in Commodites, n in Noeuds;n!=k], sum(f[k,j=>n] for j in Voisins[n]) - sum(f[k,n=>j] for j in Voisins[n]) >= 0)
            @constraint(m, source_k[k in Commodites], sum(f[k, k=>j] for j in Voisins[k]) - sum(f[k,j=>k] for j in Voisins[k]) == x[k]-r[k])

            if is_non_reserve
                @constraint(m, flux_arc_1_nr[k in Commodites, d in union(Arcs,ArcsFictif)], g[k,d] <= v[d])
                @constraint(m, flux_arc_2_nr[k in Commodites, d in union(Arcs,ArcsFictif)], g[k,d] <= 1-x[k]) # redondante ?
                @constraint(m, conservation_flux_1_nr[k in Commodites], sum(g[k,alpha=>j] for j in VoisinsFictif[alpha]) - sum(g[k,j=>alpha] for j in VoisinsFictif[alpha]) <= 1)
                @constraint(m, conservation_flux_2_nr[k in Commodites,  n in Noeuds;n!=k], sum(g[k,j=>n] for j in VoisinsFictif[n]) - sum(g[k,n=>j] for j in VoisinsFictif[n]) == 0)
                @constraint(m, sink_k_nr[k in Commodites], sum(g[k, j=>k] for j in VoisinsFictif[k]) - sum(g[k,k=>j] for j in VoisinsFictif[k]) == 1 - x[k])
            end
        end

        return m
    end
