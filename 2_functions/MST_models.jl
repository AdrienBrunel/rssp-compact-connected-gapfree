# ==============================================================================
# 0 - DEV
# ==============================================================================
    function ReserveSiteSelection_SpatialConstraints(instance, gridgraph, params, is_beta, is_non_reserve, is_callbacks, is_damier, is_rmax)

        ## LECTURE DES DONNEES -------------------------------------------------
        Noeuds               = gridgraph.Noeuds
        Voisins              = gridgraph.Voisins
        Aretes               = gridgraph.Aretes
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

        if is_damier
            NoeudsNoirs  = gridgraph.NoeudsNoirs
            NoeudsBlancs = gridgraph.NoeudsBlancs
        end

        if is_non_reserve
            NoeudsFictif   = gridgraph.NoeudsFictif
            VoisinsFictif  = gridgraph.VoisinsFictif
            AretesFictif   = gridgraph.AretesFictif
            ArcsFictif     = gridgraph.ArcsFictif
            alpha          = gridgraph.alpha
            N_noeudsfictif = length(NoeudsFictif)
        end

        if is_rmax
            Rmax = params.Rmax
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

        if is_damier
            Commodites = NoeudsNoirs
        else
            Commodites = Noeuds
        end
		if is_beta
        	@variable(m, z[Arcs], Bin);
		end
        @variable(m, u[Arcs], Bin);
        @variable(m, r[Noeuds], Bin);
        if is_non_reserve
            @variable(m, x[union(Noeuds,alpha)], Bin);
            @variable(m, v[union(Arcs,ArcsFictif)], Bin);
            @variable(m, g[Commodites,union(Arcs,ArcsFictif)], Bin);
        else
            @variable(m, x[Noeuds], Bin);
        end

        if is_rmax
            dmin = shortest_distances(N_noeuds,Voisins)
            NoeudsProche  = Vector{Vector{Int}}()
            VoisinsProche = Vector{Vector{Vector{Int}}}()
            ArcsProche    = Vector{Vector{Pair{Int,Int}}}()
            for k in Noeuds
                push!(NoeudsProche, Vector{Int}())
                push!(VoisinsProche, Vector{Vector{Int}}())
                push!(ArcsProche, Vector{Pair{Int,Int}}())
                for j in Noeuds
                    push!(VoisinsProche[k], Vector{Int}())
                    if dmin[k,j]<=Rmax
                        push!(NoeudsProche[k],j)
                        for i in Voisins[j]
                            if dmin[k,i]<=(Rmax-1)
                                push!(VoisinsProche[k][j],i)
                            end
                        end
                    end
                end
            end
            @variable(m, f[k in Commodites, i in NoeudsProche[k], j in Voisins[i];dmin[k,i]<=Rmax-1], Bin);
        else
            @variable(m, f[Commodites,Arcs], Bin);
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

		if is_beta
	        @constraint(m, linzi[d in Arcs],  z[d] - x[d[1]] <= 0)
        	@constraint(m, linzk[d in Arcs],  z[d] - x[d[2]] <= 0)
		end

        # Arbre couvrant [partie commune à tous les modèles]
        @constraint(m, arc_unique[a in Aretes],  u[a[1]=>a[2]] + u[a[2]=>a[1]] <= 1)
        @constraint(m, cycle, sum(u[d] for d in Arcs) == sum(x[j] for j in Noeuds)-1)
        @constraint(m, arc_reserve_1[d in Arcs],  u[d] <= x[d[1]])
        @constraint(m, arc_reserve_2[d in Arcs],  u[d] <= x[d[2]])
        @constraint(m, racine_reserve[j in Noeuds],  r[j] <= x[j])
        @constraint(m, racine_unique, sum(r[j] for j in Noeuds) <= 1)

        # Arbre couvrant de la non-reserve
        if is_non_reserve
            @constraint(m, locked_out, x[alpha] == 0)
            @constraint(m, arc_unique_nr[a in union(Aretes,AretesFictif)],  v[a[1]=>a[2]] + v[a[2]=>a[1]] <= 1)
            @constraint(m, cycle_nr, sum(v[d] for d in union(Arcs,ArcsFictif)) == sum(1-x[j] for j in union(Noeuds,alpha))-1)
            @constraint(m, arc_reserve_1_nr[d in union(Arcs,ArcsFictif)],  v[d] <= 1-x[d[1]])
            @constraint(m, arc_reserve_2_nr[d in union(Arcs,ArcsFictif)],  v[d] <= 1-x[d[2]])
        end

        # Visualisation : un seul arc sélectionné par voisin dans la réserve (si ce n'est pas la racine)
        @constraint(m, visu_r[j in Noeuds], sum(u[j=>i] for i in Voisins[j]) >=  x[j]-r[j])

        # Réduction damier
        if is_damier
            NoeudsBlancsInterieurs = setdiff(NoeudsBlancs,NoeudsPeripheriques)
            @constraint(m, damier_1[j in NoeudsBlancsInterieurs], 1-x[j] <=  sum(1-x[i] for i in Voisins[j]))
            @constraint(m, damier_2[j in NoeudsBlancs], x[j] <=  sum(x[i] for i in Voisins[j]))
        end

        if is_rmax
            @constraint(m, source_rmax[j in Noeuds, i in Noeuds;dmin[i,j]>Rmax], x[i] <=  1-r[j])
        end


        # Contrainte sur les variables de flux
        if is_callbacks

            function my_callback_function(cb_data)

                if is_non_reserve
                    x_val = zeros(Int,N_noeudsfictif)
                    for j in NoeudsFictif
                        val = callback_value(cb_data, x[j])
                        if abs(val-round(Int,val)) < 1e-4
                            x_val[j] = round(Int,val)
                        else
                            return
                        end
                    end
                else
                    x_val = zeros(Int,N_noeuds)
                    for j in Noeuds
                        val = callback_value(cb_data, x[j])
                        if abs(val-round(Int,val)) < 1e-4
                            x_val[j] = round(Int,val)
                        else
                            return
                        end
                    end
                end

                # Lazy cuts pour la connectivité de la réserve
                Reserve_Components = connected_components(Voisins,x_val)
                if length(Reserve_Components) > 1

                    println("$(length(Reserve_Components)) Composante connexes de la réserve : $(Reserve_Components)")

                    for c in 1:length(Reserve_Components)

                        # choix de la commodité
                        Component = intersect(Reserve_Components[c],Commodites)
                        k = Component[findmax(Rentability[Component])[2][1]]
                        println("La commodité choisie pour la réserve est $(k)")

                        # Ajout des contraintes de flux de la réserve
                        if is_rmax
                            for i in NoeudsProche[k]
                                if dmin[k,i]==Rmax
                                    continue
                                end
                                for j in Voisins[i]
                                    flux_arc_1 = @build_constraint(f[k,i,j] <= u[i=>j])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_1)
                                    flux_arc_2 = @build_constraint(f[k,i,j] <= x[k])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_2)
                                end
                            end
                            source_k = @build_constraint(sum(f[k,k,j] for j in Voisins[k]) - sum(f[k,j,k] for j in Voisins[k]) == x[k]-r[k])
                            MOI.submit(m, MOI.LazyConstraint(cb_data), source_k)

                            for i in NoeudsProche[k]
                                if i==k
                                    continue
                                end
                                if dmin[k,i] == Rmax
                                    conservation_flux_1 = @build_constraint(sum(f[k,j,i] for j in VoisinsProche[k][i]) == r[i])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1)
                                else
                                    conservation_flux_1 = @build_constraint(sum(f[k,j,i] for j in VoisinsProche[k][i]) - sum(f[k,i,j] for j in VoisinsProche[k][i]) <= r[i])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1)
                                    conservation_flux_2 = @build_constraint(sum(f[k,j,i] for j in VoisinsProche[k][i]) - sum(f[k,i,j] for j in VoisinsProche[k][i]) >= 0)
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_2)
                                end
                            end
                        else
                            for i in Noeuds
                                for j in Voisins[i]
                                    flux_arc_1 = @build_constraint(f[k,i=>j] <= u[i=>j])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_1)
                                    flux_arc_2 = @build_constraint(f[k,i=>j] <= x[k])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_2)
                                end
                            end
                            source_k = @build_constraint(sum(f[k,k=>j] for j in Voisins[k]) - sum(f[k,j=>k] for j in Voisins[k]) == x[k]-r[k])
                            MOI.submit(m, MOI.LazyConstraint(cb_data), source_k)

                            for i in Noeuds
                                if i != k
                                    conservation_flux_1 = @build_constraint(sum(f[k,j=>i] for j in Voisins[i]) - sum(f[k,i=>j] for j in Voisins[i]) <= r[i])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1)
                                    conservation_flux_2 = @build_constraint(sum(f[k,j=>i] for j in Voisins[i]) - sum(f[k,i=>j] for j in Voisins[i]) >= 0)
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_2)
                                end
                            end
                        end
                    end
                else
                    if is_non_reserve

                        # Lazy cuts pour la connectivité de la non-réserve
                        NonReserve_Components = connected_components(VoisinsFictif,1 .- x_val)
                        if length(NonReserve_Components) > 1

                            println("$(length(NonReserve_Components)) Composante Connexes de la Non-Réserve: $(NonReserve_Components)")

                            for c in 1:length(NonReserve_Components)
                                if !(sum(findall(NonReserve_Components[c] .== alpha))>1)

                                    # choix de la commodité
                                    Component = intersect(NonReserve_Components[c],Commodites)
                                    k = Component[findmin(Rentability[Component])[2][1]]
                                    println("La commodité choisie pour la non-réserve est $(k)")

                                    # Ajout des contraintes de flux de la non-réserve
                                    for d in Arcs
                                        flux_arc_1_nr = @build_constraint(g[k,d] <= v[d])
                                        MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_1_nr)
                                        flux_arc_2_nr = @build_constraint(g[k,d] <= 1-x[k])
                                        MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_2_nr)
                                    end
                                    source_k_nr = @build_constraint(sum(g[k, k=>j] for j in VoisinsFictif[k]) - sum(g[k,j=>k] for j in VoisinsFictif[k]) == 1 - x[k])
                                    MOI.submit(m, MOI.LazyConstraint(cb_data), source_k_nr)

									conservation_flux_1_nr = @build_constraint(sum(g[k,j=>alpha] for j in VoisinsFictif[alpha]) - sum(g[k,alpha=>j] for j in VoisinsFictif[alpha]) <= 1)
                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1_nr)

									for i in Noeuds
		                                if i != k
											conservation_flux_2_nr = @build_constraint(sum(g[k, i=>j] for j in VoisinsFictif[i]) - sum(g[k,j=>i] for j in VoisinsFictif[i]) == 0)
		                                    MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_2_nr)
										end
									end

									#@constraint(m, conservation_flux_2_nr[k in Noeuds,  n in Noeuds;n!=k], sum(g[k,j=>n] for j in VoisinsFictif[n]) - sum(g[k,n=>j] for j in VoisinsFictif[n]) == 0)


                                end
                            end
                        end
                    end
                end

                if is_rmax

					r_val = zeros(Int,N_noeuds)
                    for j in Noeuds
						val = callback_value(cb_data, r[j])
						if abs(val-round(Int,val)) < 1e-4
                            r_val[j] = round(Int,val)
                        else
                            return
                        end
                    end

                    Reserve = Noeuds[findall(x_val .== 1)]

					Voisins_InReserve = Vector{Vector{Int}}()
				    for k in Reserve
				        push!(Voisins_InReserve, Vector{Int}())
				    end

				    cpt = 0
				    for n1 in Reserve
				        cpt = cpt+1
				        for n2 in Voisins[n1]
				            if sum(n2 .== Reserve) > 0
				                push!(Voisins_InReserve[cpt],findall(n2 .== Reserve)[1])
				            end
				        end
				    end
				    dmin_in_reserve = shortest_distances(length(Reserve),Voisins_InReserve)
					idx_source = findall(r_val .== 1)[1]
					Reserve_TooFar = Reserve[findall(dmin_in_reserve[findall(idx_source .== Reserve)[1],:] .> Rmax)]

                    if length(Reserve_TooFar) > 0
						for k in intersect(Reserve_TooFar,Commodites)

							#if (dmin_in_reserve[findall(idx_source .== Reserve)[1],findall(k .== Reserve)[1]] == maximum(dmin_in_reserve))

								for i in NoeudsProche[k]
									if dmin[k,i]==Rmax
										continue
									end
									for j in Voisins[i]
										flux_arc_1 = @build_constraint(f[k,i,j] <= u[i=>j])
										MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_1)
										flux_arc_2 = @build_constraint(f[k,i,j] <= x[k])
										MOI.submit(m, MOI.LazyConstraint(cb_data),flux_arc_2)
									end
								end
								source_k = @build_constraint(sum(f[k,k,j] for j in Voisins[k]) - sum(f[k,j,k] for j in Voisins[k]) == x[k]-r[k])
								MOI.submit(m, MOI.LazyConstraint(cb_data), source_k)

								for i in NoeudsProche[k]
									if i==k
										continue
									end
									if dmin[k,i] == Rmax
										conservation_flux_1 = @build_constraint(sum(f[k,j,i] for j in VoisinsProche[k][i]) == r[i])
										MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1)
									else
										conservation_flux_1 = @build_constraint(sum(f[k,j,i] for j in VoisinsProche[k][i]) - sum(f[k,i,j] for j in VoisinsProche[k][i]) <= r[i])
										MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_1)
										conservation_flux_2 = @build_constraint(sum(f[k,j,i] for j in VoisinsProche[k][i]) - sum(f[k,i,j] for j in VoisinsProche[k][i]) >= 0)
										MOI.submit(m, MOI.LazyConstraint(cb_data),conservation_flux_2)
									end
								end

								rmax_inreserve = @build_constraint(sum(f[key] for key in eachindex(f) if key[1]==k) <= Rmax)
								MOI.submit(m, MOI.LazyConstraint(cb_data),rmax_inreserve)
							#end
						end
                    end
                end
            end
            MOI.set(m, MOI.LazyConstraintCallback(), my_callback_function)

        else
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
