using Pkg
#Pkg.add("EcologicalNetworks")

#import Pkg; Pkg.add("CSV"); Pkg.add("DataFrames")
using CSV
using DataFrames
using EcologicalNetworks

#cd("/home/katherine/Documents/Collaborations/Indian_foodwebs")




webnames = ["Premonsoon","Postmonsoon","Monsoon"]

function openwebs(webnames)
	dfs = []

	for web in webnames
		filename = "data/Tampara/" * web * "_adjacency_matrix.csv"

		df = CSV.read(filename,DataFrame)

		# Append the df to the list
        push!(dfs, df)

    end

    return dfs

end

webs = openwebs(webnames)

function getspecies(webs)
	Ss = []

	for web in webs
		# Get the first column which has the names of all the species
		S = web.Column1

		# Append S to the list
        push!(Ss, S)

    end

    return Ss

end


Ss = getspecies(webs)


function convertWebsToBoolArray(webs)
	As = []

	for web in webs
		## Remove the column with the names
		A = Array(web[:,2:end])

		# Convert to boolean
		A = isone.(A)

		# Append A to the list
        push!(As, A)

    end

    return As

end

As = convertWebsToBoolArray(webs)


Ns = [UnipartiteNetwork(As[j],Ss[j]) for j in 1:length(As)]

links_byweb = [links(N) for N in Ns]
C_byweb = [connectance(N) for N in Ns]


function create_dict_array(dicts, webnames, Ss)
	species = collect(keys(dicts[1]))

	for d in dicts[2:end]
	    species = union(species, keys(d))
	end
	
	#data = Array[length(species) + 1, length(dicts)]
	data = Array{Any}(undef, length(species)+1, length(dicts) + 1)

	data[1, 1] = "species"
	
	for i in 2:length(species)+1
	    data[i, 1] = species[i-1]
	end
	
	for j in 1:length(dicts)
	    for i in 2:length(species)+1
	        data[i, j+1] = get(dicts[j], species[i-1], "")
	    end
	end

	data[1,2:4] = webnames

	return data
end


degrees_byweb = [EcologicalNetworks.degree(N) for N in Ns]
degree_array = create_dict_array(degrees_byweb, webnames)
### Write to csv
#CSV.write("data/Tampara/processed/degree_dataframe.csv",  Tables.table(degree_array), writeheader=false)



outdegrees_byweb = [EcologicalNetworks.degree(N, dims = 1) for N in Ns]
outdegree_array = create_dict_array(outdegrees_byweb, webnames)
### Write to csv
#CSV.write("data/Tampara/processed/outdegree_dataframe.csv",  Tables.table(outdegree_array), writeheader=false)


specificity_byweb = [specificity(N) for N in Ns]
specificity_array = create_dict_array(specificity_byweb, webnames)
### Write to csv
#CSV.write("data/Tampara/processed/specificity_dataframe.csv",  Tables.table(specificity_array), writeheader=false)


centrality_degree_byweb = [centrality_degree(N) for N in Ns]
centrality_degree_array = create_dict_array(centrality_degree_byweb, webnames)
### Write to csv
#CSV.write("data/Tampara/processed/centrality_degree_dataframe.csv",  Tables.table(centrality_degree_array), writeheader=false)


centrality_closeness_byweb = [centrality_closeness(N) for N in Ns]
centrality_closeness_array = create_dict_array(centrality_closeness_byweb, webnames)
### Write to csv
#CSV.write("data/Tampara/processed/centrality_closeness_dataframe.csv",  Tables.table(centrality_closeness_array), writeheader=false)


overlap_byweb = [overlap(N) for N in Ns]       # calculate overlap based on prey (dims = 1) or predators (dims = 2)
#overlap_array = create_dict_array(overlap_byweb, webnames)   # This is more complicated because its done by pairs of species

AJS_byweb = [AJS(N) for N in Ns]
#AJS_array = create_dict_array(AJS_byweb, webnames) # This is more complicated because its done by pairs of species

trophic_level_byweb = [trophic_level(N) for N in Ns]
trophic_level_array = create_dict_array(trophic_level_byweb, webnames)
### Write to csv
#CSV.write("data/Tampara/processed/trophic_level_dataframe.csv",  Tables.table(trophic_level_array), writeheader=false)


omnivory_byweb = [omnivory(N) for N in Ns]
omnivory_array = create_dict_array(omnivory_byweb, webnames)
### Write to csv
#CSV.write("data/Tampara/processed/omnivory_dataframe.csv",  Tables.table(omnivory_array), writeheader=false)


## Motifs

## list all the motifs
#unipartitemotifs()

# Get motif counts

## Get all S1 motifs in the first network
#find_motif( Ns[1], unipartitemotifs().S1)

#motif_counts = [size(find_motif(N,m), 1) for m in unipartitemotifs() for N in Ns]

# Get all the motif tuples for each motif and each web
motif_lists = [find_motif(N,m) for m in unipartitemotifs(), N in Ns]

# motif 3 doesn't seem to exist in any of our webs, so this breaks if we do it normally
motifs = [1,2,4,5,6,7,8,9,10,11,12,13]

# Make a dataframe to store everything
species_motif_counts_df = DataFrame(web = String[], motif = Symbol[], species = String[], count = Int[])
# Get the motif names
motifnames = fieldnames(typeof(unipartitemotifs()))


# Run through each network, motif, and species in each network 
for N in 1:size(webnames,1)
	for m in motifs
		# Remember to use Ss[N] because some species don't occur in all networks
		for sp in 1:size(Ss[N],1)
			# get the name of the web, motif, and species
			web = webnames[N]
			motif = motifnames[m]
			species = Ss[N][sp]

			# Count how often that species appears in that motif
			# This could be edited to get each motif position, by adding a value in the final position of motif_lists. E.g. motif_lists[m,N][x][1,position]
			count = sum([sp in motif_lists[m,N][x][1,] for x in 1:size(motif_lists[m,N],1)]) 

			# Add a row to the dataframe
			species_motif_counts_df = push!(species_motif_counts_df, [web, motif, species, count])
		end
	end
end



CSV.write("data/Tampara/processed/species_motif_counts.csv", species_motif_counts_df, writeheader=FALSE)

