### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8db63e50-666f-11f1-97ea-bf31267d99db
begin
	using PlutoUI;
	PlutoUI.TableOfContents(title = "Gyroscopic Stabilization with Reaction Wheels")
end

# ╔═╡ 7b24588d-d6ba-415e-9cf1-e05fd4df864c
begin
	import LinearAlgebra, Symbolics,  ControlSystemsBase
	using Latexify
end

# ╔═╡ cef7f077-28d2-4363-8787-0b384212dc1a
md"""
# Simulation of Controllability of Gyrostat
### 1. Configuration
#### Mechanical
Rod Mass: $(@bind rodMass Slider(1:0.1:80, default=15.0, show_value=true)) g \
Probe Mass: $(@bind probeMass Slider(1:0.1:100, default=50.0, show_value = true)) g \
Probe Length: $(@bind probeLength Slider(0.1:0.1:3, default=1.0, show_value = true)) cm \
Rod Length: $(@bind rodLength Slider(1:0.1:30, default=15.0, show_value = true)) cm \
Rod Diameter: $(@bind rodDiameter Slider(0.1:0.1:2, default=2.0, show_value = true)) cm \
Wheel Inertia: $(@bind wheelInertia NumberField(0.1:0.1:50, default=25.0)) g cm^2\
Spinning damping: $(@bind dampingSpin NumberField(0.1:0.1:50, default=25.0)) g cm\
Transverse damping: $(@bind dampingTransverse NumberField(0.1:0.1:50, default=25.0)) g cm\
"""

# ╔═╡ 77d79145-70db-44e2-9f20-c0fd2e383ec7
md"""
Probe Inertia: $(probeInertia = (1/6*probeMass*probeLength^2 + probeMass*rodLength^2)) g cm^2 \
Axial Rod Inertia: $(axialInertia = 1/8 * rodMass * (rodDiameter)^2 + 1/6*probeMass*probeLength^2 + wheelInertia) g cm^2 \
Transverse Rod Inertia: $(transverseInertia = 1/3*rodMass*rodLength^2 + probeInertia + wheelInertia) g cm^2\
Center of Mass: $(rcm = (rodMass/2 + probeMass)*rodLength/(rodMass + probeMass)) cm
"""

# ╔═╡ 179b51a6-05f8-4810-b68c-b330ff1de0e1
begin
	physics = (
		I_T    = transverseInertia * 1e-7,   # g*cm^2 to kgm^2
	    I_A    = axialInertia * 1e-7,     # g*cm^2 to kgm^2
	    J_s    = wheelInertia * 1e-7,      # g*cm^2 to kg*m^2
		M_total = (probeMass + rodMass)/1000,     # g to kg
		g 		= 9.810,    # m/s^2
		r_cm = rcm/100,
		b_perp = dampingTransverse*1e-5,
		b_axial = dampingSpin*1e-5,
	)
end

# ╔═╡ affe95c5-4116-445b-9cd5-4b4f179e8d97
begin
	function GenerateStateSpace(phys, sens)
		I_A = phys.I_A
		I_T = phys.I_T
		b_C = phys.b_axial
		b_T = phys.b_perp
		grav = phys.g
		mass = phys.M_total
		J_s = phys.J_s
		rcm = phys.r_cm
		ω = [sens[1], sens[2], sens[3]]
		Γ = [sens[4], sens[5], sens[6]]
		u = [sens[7], sens[8], sens[9]]

		#develop tangent basis vectors for plane
	
		a = abs(Γ[1]) < 0.95 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0] # choose v1 such that it can never be zero after gram schmidt.
		v1 = LinearAlgebra.normalize(a - LinearAlgebra.dot(a, Γ))*Γ
		v2 = LinearAlgebra.cross(Γ, v1) #v3 cross v1 equals v2 in rhr orthonormal
		T = hcat(v1, v2)
		
		Im = [I_T 0 0; 0 I_T 0; 0 0 I_A]
		A_ωω = Im\([-b_T -I_A*ω[3]+I_T*ω[3]-J_s*u[3] -I_A*ω[2]+I_T*ω[2]+J_s*u[2];
			   I_A*ω[3]-I_T*ω[3]+J_s*u[3] -b_T I_A*ω[1]-I_T*ω[1]-J_s*u[1];
			   -J_s*u[2] J_s*u[1] -b_C])
		A_ωΓ = [0 mass*grav*rcm/I_T 0;
			   -mass*grav*rcm/I_T 0 0;
			   0 0 0]*T
		A_ωu = Im\(J_s*S(ω))
		
		A_Γω = T'*S(Γ)
		A_ΓΓ = T'*(-S(ω))*T
		A_Γu = T'*[0 0 0;
				   0 0 0;
				   0 0 0]
		
		A_uω = LinearAlgebra.zeros(3,3)
		A_uΓ = LinearAlgebra.zeros(3,3)*T
		A_uu = LinearAlgebra.zeros(3,3)
		
		A = [A_ωω A_ωΓ A_ωu;
			A_Γω A_ΓΓ A_Γu;
			A_uω A_uΓ A_uu]
		
		B = [-J_s/I_T 0 0;
			  0 -J_s/I_T 0;
			  0 0 -J_s/I_A;  # sign fix
			  0 0 0; 
			  0 0 0;  
			  1 0 0; 
			  0 1 0; 
			  0 0 1]

		return A, B
	end
end

# ╔═╡ 40d46ad9-23cf-4e53-a37b-d06615dd49a0
begin
	function GenerateReferenceVector(phys, sens)
		#using euler coordinates theta (nutation), phi (precession), gamma (spin)
		#solve the steady precession, steady spin, at fixed nutation angle.
		#with no damping, the angular momentum of the system would keep the rod at fixed nutation with gravity unable to overcome, meaning the transverse inputs can be set to zero and the coincident input gamma with rod gamma can be used for stabilization.

		#the solution will yield a quadratic in phi, whose determinant describes the necessary velocity u_gamma and spin need to be at.

		I_A = phys.I_A
		I_T = phys.I_T
		b_C = phys.b_axial
		b_T = phys.b_perp
		grav = phys.g
		mass = phys.M_total
		J_s = phys.J_s
		rcm = phys.r_cm
		
		ω = [sens[1], sens[2], sens[3]]
		Γ = [sens[4], sens[5], sens[6]]
		u = [sens[7], sens[8], sens[9]]
		
		#gamma axis angular momentum
		H_γ = I_A*
	end
end

# ╔═╡ bf736020-ddce-4045-bc8c-2072d4aad86c
begin
	function Control_lqr(phys, sens)
		
		Χ_ref = GenerateReferenceVector(phys, sens)
		A, B = GenerateStateSpace(phys, sens)
		K = ControlSystemsBase.lqr(A, B, Q, R)
		return output_i, output_j, output_γ
	end
	p_lqr = (
		Physics = physics,
		ControlLaw = (phys, sensor_vec) -> Control_lqr(phys, sensor_vec)
	)
end

# ╔═╡ 7ec6f4d3-a5bf-4cab-8fbf-c8f2034c3b69
begin
	function Control_adrc(phys, ref, sens)
		A, B = GenerateStateSpace(phys, sens)
		return output_i, output_j, output_γ
	end
	p_adrc = (
		Physics = physics,
		ControlLaw = (phys, sensor_vec) -> Control_adrc(phys, sensor_vec)
	)
end

# ╔═╡ e29f87fb-503a-4bf7-ba49-d37950c2602b
begin
	function GyrostatPhysics!(du, u, p, t)
		#populate data
		ω_i, ω_j, ω_γ, q_0, q_1, q_2, q_3, Ω_i, Ω_j, Ω_γ = u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8], u[9], u[10]
		
		#simulate sensor reads
		Γ_i = 2*(q1*q3 - q0*q2)
		Γ_j = 2*(q2*q3 + q0*q1)
		Γ_γ = q_0^2+q_3^2-q_2^2-q_1^2
		
		sensor = [ω_i, ω_j, ω_γ, Γ_i, Γ_j, Ω_i, Ω_j, Ω_γ]
		
		#control law
		dΩ_i, dΩ_j, dΩ_γ = p.ControlLaw(p.Physics, sensor)
		
		#update	
		MassMatrix = [p.Physics.I_T 0 0; 
					  0 p.Physics.I_T 0; 
					  0 0 p.Physics.I_A]
		GravityVector = 2*p.Physics.M_total*p.Physics.g*p.Physics.r_cm*[q0*q1+q2*q3;
												q0*q2-q1*q3;
												0]
		Beta = [p.Physics.b_perp; 
				p.Physics.b_perp; 
				p.Physics.b_axial]
		
		Skew3 = [0 -ω_γ ω_j; 
				 ω_γ 0 -ω_i; 
				 -ω_j ω_i 0]
		Skew4 = [0 -ω_i -ω_j -ω_γ; 
				 ω_i 0 ω_γ -ω_j; 
				 ω_j -ω_γ 0 ω_i; 
				 ω_γ ω_j -ω_i 0]
		
		du[4] = (1/2)*(Skew4[1,1]*q0 + Skew4[1,2]*q1 + Skew4[1,3]*q2 + Skew4[1,4]*q3)
		du[5] = (1/2)*(Skew4[2,1]*q0 + Skew4[2,2]*q1 + Skew4[2,3]*q2 + Skew4[2,4]*q3)
		du[6] = (1/2)*(Skew4[3,1]*q0 + Skew4[3,2]*q1 + Skew4[3,3]*q2 + Skew4[3,4]*q3)
		du[7] = (1/2)*(Skew4[4,1]*q0 + Skew4[4,2]*q1 + Skew4[4,3]*q2 + Skew4[4,4]*q3)

		rhs = (GravityVector .- Beta .* ω_vec .- Skew3 * (MassMatrix * [ω_i; ω_j; ω_γ] .+ p.Physics.J_s.*[ω_i; ω_j; ω_γ]) .- p.Physics.J_s.*[dΩ_i; dΩ_j; dΩ_γ])
		
		du[1] = rhs[1] / MassMatrix[1,1]
		du[2] = rhs[2] / MassMatrix[2,2]
		du[3] = rhs[3] / MassMatrix[3,3]

		du[8] = dΩ_i
		du[9] = dΩ_j
		du[10] = dΩ_γ
	end
end

# ╔═╡ 768d624b-44a3-4053-8318-4c2ab32a0f7d
begin 

end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ControlSystemsBase = "aaaaaaaa-a6ca-5380-bf3e-84a91bcd477e"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
ControlSystemsBase = "~1.20.3"
Latexify = "~0.16.10"
PlutoUI = "~0.7.80"
Symbolics = "~7.19.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "2da7d161b5a596c37d0eb9a65bcf59a78b8139f9"

[[deps.ADTypes]]
git-tree-sha1 = "bbc22a9a08a0ef6460041086d8a7b27940ed4ffd"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.22.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "2eeb2c9bef11013efc6f8f97f32ee59b146b09fb"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.44"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "0761717147821d696c9470a7a86364b2fbd22fd8"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.2"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "54f895554d05c83e3dd59f6a396671dae8999573"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.24.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.ControlSystemsBase]]
deps = ["ForwardDiff", "Hungarian", "LinearAlgebra", "MacroTools", "MatrixEquations", "MatrixPencils", "Polynomials", "PrecompileTools", "Printf", "Random", "RecipesBase", "StaticArraysCore", "UUIDs"]
git-tree-sha1 = "712db308e3ad3e974bd2df20d53732554270c994"
uuid = "aaaaaaaa-a6ca-5380-bf3e-84a91bcd477e"
version = "1.20.3"

    [deps.ControlSystemsBase.extensions]
    ControlSystemsBaseDSPExt = ["DSP"]
    ControlSystemsBaseImplicitDifferentiationExt = ["ImplicitDifferentiation", "ComponentArrays"]
    ControlSystemsBaseMakieExt = ["Makie"]

    [deps.ControlSystemsBase.weakdeps]
    ComponentArrays = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
    DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
    ImplicitDifferentiation = "57b37032-215b-411a-8a7c-41a003a55207"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "a55766a9c8f66cf19ffcdbdb1444e249bb4ace33"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.4.6"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "FunctionMaps", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4599e0cd684f3ff6cbbab73c77553a3d01a8d74d"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.18"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "StarAlgebras", "Test"]
git-tree-sha1 = "5bfabc3827dfdd164359bad6800c115a81280c00"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.6"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "afb7c51ac63e40708a3071f80f5e84a752299d4f"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.39"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FunctionMaps]]
deps = ["CompositeTypes", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "31bd99a57edf98990d1c21486032963955450e8d"
uuid = "a85aefff-f8ca-4649-a888-c8e5398bc76c"
version = "0.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7eb45fe833a5b7c51cf6d89c5a841d5967e44be3"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.14.0"

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

    [deps.Graphs.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
    SharedArrays = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Hungarian]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4f84db415ccb0ea750b10738bfecdd55388fd1b6"
uuid = "e91730f6-4275-51fb-a7a0-7064cfbd3b39"
version = "0.7.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LinearMaps]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7f6be2e4cdaaf558623d93113d6ddade7b916209"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.4"

    [deps.LinearMaps.extensions]
    LinearMapsChainRulesCoreExt = "ChainRulesCore"
    LinearMapsSparseArraysExt = "SparseArrays"
    LinearMapsStatisticsExt = "Statistics"

    [deps.LinearMaps.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixEquations]]
deps = ["LinearAlgebra", "LinearMaps"]
git-tree-sha1 = "51f3fade0b4ff2cf90b36b3312425460631abb56"
uuid = "99c1a7ee-ab34-5fd5-8076-27c950a045f4"
version = "2.5.6"

[[deps.MatrixPencils]]
deps = ["LinearAlgebra", "MatrixEquations", "Polynomials", "Random"]
git-tree-sha1 = "d852fb5efe6eea8f764d44a3e6a559664fa88e20"
uuid = "48965c70-4690-11ea-1f13-43a2532b2fa8"
version = "1.9.1"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics", "StarAlgebras"]
git-tree-sha1 = "4838893d9b035c2f6967c0d533350e1755b58a70"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.19"

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

    [deps.MultivariatePolynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "7c25249fc13a070f5ba433c50e21e22bb33c6fb0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.7.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "fbc875044d82c113a9dee6fc14e16cf01fd48872"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.80"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "2d99b4c8a7845ab1342921733fa29366dae28b24"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"
    PolynomialsRecipesBaseExt = "RecipesBase"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cfcdc949c4660544ab0fdeed169561cb22f835f4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.18"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2700b235561b0335d5bef7097a111dc513b8655e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.2"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "MutableArithmetics", "SparseArrays"]
git-tree-sha1 = "235b1f9d287bbf34083b3d0829343a7942c0ad1c"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.3.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "5085671d2cba1eb02136a3d6661c583e801984c1"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.1.0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "Dictionaries", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "Graphs", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "98c23082cba20c6ae17e920600c89e4f6491e618"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.24.2"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsChainRulesCoreExt = "ChainRulesCore"
    SymbolicUtilsDistributionsExt = "Distributions"
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "c0ff31ee7be485d05d15a0b32a2cd615e90cbfc0"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.19.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsHypergeometricFunctionsExt = "HypergeometricFunctions"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    HypergeometricFunctions = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═8db63e50-666f-11f1-97ea-bf31267d99db
# ╠═cef7f077-28d2-4363-8787-0b384212dc1a
# ╠═77d79145-70db-44e2-9f20-c0fd2e383ec7
# ╠═7b24588d-d6ba-415e-9cf1-e05fd4df864c
# ╠═179b51a6-05f8-4810-b68c-b330ff1de0e1
# ╟─affe95c5-4116-445b-9cd5-4b4f179e8d97
# ╠═40d46ad9-23cf-4e53-a37b-d06615dd49a0
# ╠═bf736020-ddce-4045-bc8c-2072d4aad86c
# ╠═7ec6f4d3-a5bf-4cab-8fbf-c8f2034c3b69
# ╠═e29f87fb-503a-4bf7-ba49-d37950c2602b
# ╠═768d624b-44a3-4053-8318-4c2ab32a0f7d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
