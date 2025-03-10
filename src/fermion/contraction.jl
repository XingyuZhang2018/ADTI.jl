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

function next_neighbor1_contraction(h1, h2, A11, A11dag, A12, A12dag, A21, A21dag, A22, A22dag, SdD, SDD,
	ACu, ARu, FLu, FLo, FRu, FRo, ACd, ARd)
	@tensoropt e = SDD[3 2; 1 10] * A11[4 3 11 12 5] * SdD[24 21; 11 12] * SdD[17 18; 20 19] * 
				   SDD[19 22; 23 21] * SdD[20 29; 30 22] * SdD[31 23; 24 25] * A11dag[10 16 17 28 18] * h1[31 30; 71] * 
				   SDD[7 6; 5 13] * A12[8 7 14 15 9] * SDD[26 35; 27 15] * A12dag[13 25 14 34 26] * 
				   SDD[40 28; 38 39] * A21[29 40 41 54 42] * SDD[53 55; 56 54] * A21dag[39 51 41 52 53] * 
				   SdD[32 34; 37 36] * SdD[47 43; 33 42] * SDD[44 36; 43 48] * SdD[59 48; 47 57] * h2[71; 32 33] *
				   SdD[37 45; 49 44] * A22[35 45 49 50 46] * SDD[60 61; 62 50] * A22dag[57 56 59 58 60] * 
				   ACu[63 2 4; 64] * ARu[64 6 8; 65] * FLu[70 16 1; 63] * FLo[69 51 38; 70] * 
				   FRu[65 27 9; 66] * FRo[66 62 46; 67] * ACd[68; 69 52 55] * ARd[67; 68 58 61]

	return e
end

function next_neighbor2_contraction(h1, h2, A11, A11dag, A12, A12dag, A21, A21dag, A22, A22dag, SdD, SDD,
	ACu, ARu, FLu, FLo, FRu, FRo, ACd, ARd)
	@tensoropt e = SDD[2 25; 1 26] * A11[30 2 48 31 3] * SDD[7 32; 8 31] * A11dag[26 6 48 27 7] *
				   SDD[4 36; 3 37] * A12[43 4 49 44 5] * SDD[11 45; 12 44] * A12dag[39 10 50 40 11] *
				   SdD[58 37; 49 38] * SdD[50 38; 59 39] * SdD[55 8; 58 9] * SdD[59 9; 57 10] * h2[71; 55 57] * 
				   SDD[14 27; 13 28] * A21[32 14 51 62 15] * SdD[60 33; 51 62] * SdD[54 16; 60 15] *
				   SdD[52 34; 61 33] * SdD[61 17; 56 16] * SDD[21 35; 22 34] * A21dag[28 20 52 29 21] * h1[54 56; 71] * 
				   SDD[18 40; 17 41] * A22[45 18 53 46 19] * SDD[23 47; 24 46] * A22dag[41 22 53 42 23] *
				   ACu[63 25 30; 64] * ARu[64 36 43; 65] * FLu[70 6 1; 63] * FLo[69 20 13; 70] *
				   FRu[65 12 5; 66] * FRo[66 24 19; 67] * ACd[68; 69 29 35] * ARd[67; 68 42 47]

	return e
end

function four_site_contraction(M11, M12, M21, M22, ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo)
	@tensoropt n = FLu[9 6; 1] * ACu[1 2; 3] * M11[6 10; 7 2] * ARu[3 4; 5] * FRu[5 8; 12] * M12[7 11; 8 4] * FLo[16 13; 9] * ACd[18; 16 17] * M21[13 17; 14 10] * ARd[20; 18 19] * FRo[12 15; 20] * M22[14 19; 15 11]
	return n
end
