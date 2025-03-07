function horizontal_contraction(h1, h2, A1, A1dag, A2, A2dag, SdD, SDD, FLo, ACu, ARu, FRo, ARd, ACd)
 	@tensoropt eh = FLo[38 7 29; 39] * ACd[37; 38 9 40] * A1dag[6 7 8 9 10] * SDD[10 40; 17 27] * SdD[8 27; 23 22] * ACu[39 30 1; 34] * A1[1 2 3 4 5] * SDD[2 30; 29 6] * SdD[21 22; 3 4] * h1[21 23; 41] * ARd[36; 37 19 32] * A2dag[16 17 18 19 20] * SDD[32 20; 14 31] * SdD[18 26; 25 16] * ARu[34 33 11; 35] * A2[11 12 13 14 15] * SDD[12 33; 5 28] * SdD[24 28; 13 26] * h2[41; 24 25] * FRo[35 31 15; 36]
	
	return eh
end

function vertical_contraction(h1, h2, A1, A1dag, A2, A2dag, SdD, SDD, ACu, FRu, FRo, ACd, FLo, FLu)
	@tensoropt ev = ACu[35 12 1; 36] * FRu[36 18 5; 37] * A1[1 2 3 4 5] * SDD[12 2; 6 11] * SdD[15 13; 3 16] * FLu[40 7 11; 35] * A1dag[6 7 8 9 10] * SDD[9 33; 23 34] * SdD[8 10; 14 13] * h1[15 14; 41] * FLo[39 24 34; 40] * A2dag[23 24 25 26 27] * SDD[27 28; 29 21] * SdD[25 31; 32 33] * ACd[38; 39 26 28] * A2[17 19 20 21 22] * SDD[17 16; 4 18] * SdD[30 19; 20 31] * h2[41; 30 32] * FRo[37 29 22; 38]

	return ev
end

# function generate_onsite_rules(;d=4,D=2,χ=20)
#     eincode = EinCode(((1,2,3,4,5),# T1

# 	(10,7,6,8,9),#T2 (dag)

#     (14,4,9,16), #swapgate(nl,nu)
#     (18,10,2,12), #swapgate(nl,nu)

#     (11,7,18,17), #E1 FLo
#     (11,12,1,13), #E2 ACu
#     (13,14,5,15), #E4 FRo
#     (17,8,16,15), #E6 ACd
# 	),
# 	(3,6) #hamiltonian (ij di dj)
# 	)

# 	size_dict = [D for _ = 1:18]
# 	size_dict[[3;6]] .= d
# 	size_dict[[11;13;15;17]] .= χ
# 	sd = Dict(i=>size_dict[i] for i = 1:18)
# 	# for seed = 1:100
# 	seed = 4
# 	Random.seed!(seed)
# 	optcode = optimize_code(eincode, sd, TreeSA())
# 	print("On site Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n")
# 	# end

# 	return optcode
# end

# function generate_next_neighbor1_rules(;d=4,D=2,χ=20)
# 	eincode = EinCode((
# 	(1,10,3,2),#SDD
# 	(4,3,11,12,5),#T11

# 	(11,21,24,12),#SdD
# 	(17,18,20,19),#SdD
# 	(19,22,23,21),#SDD
# 	(20,29,30,22),#SdD
# 	(24,23,31,25),#SdD
	
# 	(10,16,17,28,18),#T11'

# 	(5,13,7,6),#SDD
# 	(8,7,14,15,9),#T12
# 	(26,35,27,15),#SDD
# 	(13,25,14,34,26),#T12'

# 	(38,39,40,28),#SDD
# 	(29,40,41,54,42),#T21
# 	(53,55,56,54),#SDD
# 	(39,51,41,52,53),#T21'

# 	(32,36,37,34),#SdD
# 	(33,42,47,43),#SdD
# 	(43,48,44,36),#SDD
# 	(47,57,59,48),#SdD
# 	(37,44,49,45),#SdD

# 	(35,45,49,50,46),#T22
# 	(60,61,62,50),#SDD
# 	(57,56,59,58,60),#T22'

# 	(63,2,4,64),#ACu
# 	(64,6,8,65),#ARu
# 	(63,16,1,70),#FLu
# 	(70,51,38,69),#FLo
# 	(65,27,9,66),#FRu
# 	(66,62,46,67),#FRo
# 	(69,52,55,68),#ACd
# 	(68,58,61,67),#ARd
# 	),
# 	(31,32,30,33) #hamiltonian (ij di dj)
# 	)
		
# 	size_dict = [D for _ = 1:70]
# 	size_dict[[11;24;31;17;20;30;49;37;32;59;47;33]] .= d
# 	size_dict[63:70] .= χ
# 	sd = Dict(i=>size_dict[i] for i = 1:70)

# 	# for seed =40:100
# 	seed = 100
# 	Random.seed!(seed)
# 	optcode = optimize_code(eincode, sd, TreeSA())


# 	print("NN1 Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n") 
# 	# You would better try some times to make it optimal (By the for-end iteration...)
# 	# end
# 	return optcode
# end
# oc_NN1_fermion = generate_next_neighbor1_rules()

# function generate_next_neighbor2_rules(;d=4,D=2,χ=20)
# 	eincode = EinCode((
# 	(1,26,2,25),#SDD
# 	(30,2,48,31,3),#T11
# 	(7,32,8,31),#SDD
# 	(26,6,48,27,7),#T11'

# 	(3,37,4,36),#SDD
# 	(43,4,49,44,5),#T12
# 	(11,45,12,44),#SDD

# 	(39,10,50,40,11),#T12'
# 	(58,38,49,37),#SdD
# 	(59,39,50,38),#SdD
# 	(58,8,55,9),#SdD
# 	(59,9,57,10),#SdD
	
#     (13,28,14,27),#SDD

# 	(32,14,51,62,15),#T21
# 	(51,33,60,62),#SdD
# 	(54,15,60,16),#SdD
# 	(52,34,61,33),#SdD
# 	(56,16,61,17),#SdD
	
# 	(21,35,22,34),#SDD
# 	(28,20,52,29,21),#T21'

# 	(17,41,18,40),#SDD
# 	(45,18,53,46,19),#T22
# 	(23,47,24,46),#SDD
# 	(41,22,53,42,23),#T22'

# 	(63,25,30,64),#ACu
# 	(64,36,43,65),#ARu
# 	(63,6,1,70),#FLu
# 	(70,20,13,69),#FLo
# 	(65,12,5,66),#FRu
# 	(66,24,19,67),#FRo
# 	(69,29,35,68),#ACd
# 	(68,42,47,67),#ARd
# 	),
# 	(54,55,56,57) #hamiltonian (ij di dj)
# 	)

# 	size_dict = [D for _ = 1:70]
# 	size_dict[[54;60;51;55;58;49;56;61;52;57;59;50]] .= d
# 	size_dict[63:70] .= χ
# 	sd = Dict(i=>size_dict[i] for i = 1:70)

# 	# for seed =40:100
# 	seed = 100
# 	Random.seed!(seed)
# 	optcode = optimize_code(eincode, sd, TreeSA())


# 	print("NN2 Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n") 
# 	# You would better try some times to make it optimal (By the for-end iteration...)
# 	# end
# 	return optcode
# end
# oc_NN2_fermion = generate_next_neighbor2_rules()