function generate_horizontal_rules(;d=4,D=2,χ=20)
	eincode = EinCode(((1,2,3,4,5),# T1
    (3,4,21,22),#swapgate(nf,nu)

	(6,7,8,9,10),#T2 (dag)

    (23,22,8,27),#swapgate(nf,nu)
    (17,27,10,40), #swapgate(nl,nu)
    (29,6,2,30), #swapgate(nl,nu)

    (16,17,18,19,20),# T4(dag)

    (25,16,18,26),#swapgate(nf,nu)
    (13,26,24,28),#swapgate(nf,nu)
    (5,28,12,33),#swapgate(nl,nu)
    (11,12,13,14,15),# T3

    (31,14,20,32),#swapgate(nl,nu)

    (39,7,29,38), #E1 FLo
    (39,30,1,34), #E2 ACu
    (34,33,11,35), #E3 ARu
    (35,31,15,36), #E4 FRo
    (37,19,32,36), #E5 ARd
    (38,9,40,37), #E6 ACd
	),
	(21,24,23,25) #hamiltonian (ij di dj)
	)

    size_dict = [D for i = 1:40]
	size_dict[[3;8;13;18;21;23;24;25]] .= d
	size_dict[34:39] .= χ
	sd = Dict(i=>size_dict[i] for i = 1:40)
	
	# for seed = 1:100
	seed = 100
	Random.seed!(seed)
	optcode = optimize_code(eincode, sd, TreeSA())
	print("Horizontal Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n")
	# end

	return optcode
end
oc_H_fermion = generate_horizontal_rules()

function generate_vertical_rules(;d=4,D=2,χ=20)
	eincode = EinCode(((1,2,3,4,5),# T1
	(6,7,8,9,10),#T2 (dag)

	(11,6,2,12), #swapgate(nl,nu)
	(14,13,8,10), #swapgate(nf,nu)
	(3,16,15,13), #swapgate(nf,nr)
	(4,18,17,16), #swapgate(nl,nu)

	(17,19,20,21,22), #T4
	(23,24,25,26,27), #T3 (dag)

	(21,29,28,27),#swapgate(nl,nu)
	(20,31,30,19),#swapgate(nf,nr)
	(32,33,25,31),#swapgate(nf,nr)
	(23,34,9,33), #swapgate(nl,nu)

	(35,12,1,36), # ACu: E3
	(36,18,5,37), # FRu: E8
	(37,29,22,38), # FRo: E4
	(39,26,28,38),# ACd: E6
	(40,24,34,39), # FLo: E1
	(35,7,11,40) # FLu: E7
	),
	(15,30,14,32) #hamiltonian (ij di dj)
	)
		
	size_dict = [D for i = 1:40]
	size_dict[[3;8;14;15;30;32;20;25]] .= d
	size_dict[35:40] .= χ
	sd = Dict(i=>size_dict[i] for i = 1:40)

	# for seed =40:100
	seed = 100
	Random.seed!(seed)
	optcode = optimize_code(eincode, sd, TreeSA())


	print("Vertical Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n") 
	# You would better try some times to make it optimal (By the for-end iteration...)
	# end
	return optcode
end
oc_V_fermion = generate_vertical_rules()

function generate_onsite_rules(;d=4,D=2,χ=20)
    eincode = EinCode(((1,2,3,4,5),# T1

	(10,7,6,8,9),#T2 (dag)

    (14,4,9,16), #swapgate(nl,nu)
    (18,10,2,12), #swapgate(nl,nu)

    (11,7,18,17), #E1 FLo
    (11,12,1,13), #E2 ACu
    (13,14,5,15), #E4 FRo
    (17,8,16,15), #E6 ACd
	),
	(3,6) #hamiltonian (ij di dj)
	)

	size_dict = [D for _ = 1:18]
	size_dict[[3;6]] .= d
	size_dict[[11;13;15;17]] .= χ
	sd = Dict(i=>size_dict[i] for i = 1:18)
	# for seed = 1:100
	seed = 4
	Random.seed!(seed)
	optcode = optimize_code(eincode, sd, TreeSA())
	print("On site Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n")
	# end

	return optcode
end

function generate_next_neighbor1_rules(;d=4,D=2,χ=20)
	eincode = EinCode((
	(1,10,3,2),#SDD
	(4,3,11,12,5),#T11

	(11,21,24,12),#SdD
	(17,18,20,19),#SdD
	(19,22,23,21),#SDD
	(20,29,30,22),#SdD
	(24,23,31,25),#SdD
	
	(10,16,17,28,18),#T11'

	(5,13,7,6),#SDD
	(8,7,14,15,9),#T12
	(26,35,27,15),#SDD
	(13,25,14,34,26),#T12'

	(38,39,40,28),#SDD
	(29,40,41,54,42),#T21
	(53,55,56,54),#SDD
	(39,51,41,52,53),#T21'

	(32,36,37,34),#SdD
	(33,42,47,43),#SdD
	(43,48,44,36),#SDD
	(47,57,59,48),#SdD
	(37,44,49,45),#SdD

	(35,45,49,50,46),#T22
	(60,61,62,50),#SDD
	(57,56,59,58,60),#T22'

	(63,2,4,64),#ACu
	(64,6,8,65),#ARu
	(63,16,1,70),#FLu
	(70,51,38,69),#FLo
	(65,27,9,66),#FRu
	(66,62,46,67),#FRo
	(69,52,55,68),#ACd
	(68,58,61,67),#ARd
	),
	(31,32,30,33) #hamiltonian (ij di dj)
	)
		
	size_dict = [D for _ = 1:70]
	size_dict[[11;24;31;17;20;30;49;37;32;59;47;33]] .= d
	size_dict[63:70] .= χ
	sd = Dict(i=>size_dict[i] for i = 1:70)

	# for seed =40:100
	seed = 100
	Random.seed!(seed)
	optcode = optimize_code(eincode, sd, TreeSA())


	print("NN1 Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n") 
	# You would better try some times to make it optimal (By the for-end iteration...)
	# end
	return optcode
end
oc_NN1_fermion = generate_next_neighbor1_rules()

function generate_next_neighbor2_rules(;d=4,D=2,χ=20)
	eincode = EinCode((
	(1,26,2,25),#SDD
	(30,2,48,31,3),#T11
	(7,32,8,31),#SDD
	(26,6,48,27,7),#T11'

	(3,37,4,36),#SDD
	(43,4,49,44,5),#T12
	(11,45,12,44),#SDD

	(39,10,50,40,11),#T12'
	(58,38,49,37),#SdD
	(59,39,50,38),#SdD
	(58,8,55,9),#SdD
	(59,9,57,10),#SdD
	
    (13,28,14,27),#SDD

	(32,14,51,62,15),#T21
	(51,33,60,62),#SdD
	(54,15,60,16),#SdD
	(52,34,61,33),#SdD
	(56,16,61,17),#SdD
	
	(21,35,22,34),#SDD
	(28,20,52,29,21),#T21'

	(17,41,18,40),#SDD
	(45,18,53,46,19),#T22
	(23,47,24,46),#SDD
	(41,22,53,42,23),#T22'

	(63,25,30,64),#ACu
	(64,36,43,65),#ARu
	(63,6,1,70),#FLu
	(70,20,13,69),#FLo
	(65,12,5,66),#FRu
	(66,24,19,67),#FRo
	(69,29,35,68),#ACd
	(68,42,47,67),#ARd
	),
	(54,55,56,57) #hamiltonian (ij di dj)
	)

	size_dict = [D for _ = 1:70]
	size_dict[[54;60;51;55;58;49;56;61;52;57;59;50]] .= d
	size_dict[63:70] .= χ
	sd = Dict(i=>size_dict[i] for i = 1:70)

	# for seed =40:100
	seed = 100
	Random.seed!(seed)
	optcode = optimize_code(eincode, sd, TreeSA())


	print("NN2 Contraction Complexity(seed=$(seed))",OMEinsum.timespace_complexity(optcode,sd),"\n") 
	# You would better try some times to make it optimal (By the for-end iteration...)
	# end
	return optcode
end
oc_NN2_fermion = generate_next_neighbor2_rules()