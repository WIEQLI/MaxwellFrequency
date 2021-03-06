using MaxwellUtils

if hasJOcTree

#---------------------------------------------------------------------------------------------------
# prepare the param from "scratch"

function getMaxwellFreqParam(x0::Array{Float64,1},
							 n::Array{Int64,1},
							 h::Array{Float64,1},
							 Srcs::Array{Array{Float64},1},
							 Recs::Array{Array{Float64},1},
							 freq::Float64,
							 nf=2,
							 nc=2,
							 Box::Array{Float64}=[2.0  2; 2  2; 2  4],
							 linSolParam::AbstractSolver=getMUMPSsolver([],1,0),
							 doFV::Bool=true,doSE=true,ncells=[1e4;2e4;])

	S = getOTMeshFromTxRx(x0,n,h,Srcs,Recs,nf,nc,Box)
	if ( ncells[1] > nnz(S) ) || (nnz(S) > ncells[2])
		i=1
		while (i<10)
			if nnz(S) < ncells[1] # not enough cells in OcTree mesh
				nf += 1
				nc += 1
			elseif nnz(S) > ncells[2] # too many cells
				if min(nf,nc)==0
				   error("Cannot construct OcTree with number of cells between $(ncells[1]) and  $(ncells[2])")
				end
				nf = max(0,nf-1)
				nc = max(0,nc-1)
			else
				break
			end
			S = getOTMeshFromTxRx(x0,n,h,Srcs,Recs,nf,nc,Box)
			i+=1
		end
	end

	M  = getOcTreeMeshFV(S,h;x0=x0)

	nEx,nEy,nEz = getEdgeNumbering(S)

   # transmitters
	s   = zeros(Complex128,sum(M.ne),length(Srcs))
	for k=1:length(Srcs)
		s[:,k] = sparse(getEdgeIntegralOfPolygonalChain(M,Srcs[k],normalize=false))
	end

   # receivers
	obs = spzeros(sum(M.ne),length(Recs))
	for k=1:length(Recs)
		obs[:,k] = sparse(getEdgeIntegralOfPolygonalChain(M,Recs[k],normalize=true))
	end
	if doSE
		pFor = getMaxwellFreqParam(M,s,obs,[],freq,linSolParam,sensitivityMethod=:Explicit)
	else
		pFor =  getMaxwellFreqParam(M,s,obs,[],freq,linSolParam)
	end
	return pFor
end

function getMaxwellFreqParam(x0::Array{Float64,1},
						     n::Array{Int64,1},
						     h::Array{Float64,1},
						     Srcs::Array{Array{Float64},1},	# Sources+Freq
							 Recs::Array{Array{Float64},1},
							 freq::Array,
							 SrcRecMap::SparseMatrixCSC;                # maps receivers to sources
							 ProbSrcMap::SparseMatrixCSC=speye(length(Srcs)),  # maps sources to pFor's, default: single source
							 nf=2,
							 nc=2,
							 Box=[2.0  2; 2  2; 2  4],
							 linSolParam::AbstractSolver=getMUMPSsolver([],1,0),
							 doFV=true,doSE=true,
							 workerList=workers(),
							 ncells=[1e4;2e4;],  # min and max number of cells per mesh
					    	 maxProbs=ceil(Int64,(2+size(ProbSrcMap,1)/length(workerList))*ones(maximum(workerList))) # maximum number of problems per worker
							)

    workerList = intersect(workers(),workerList)
    if isempty(workerList)
    	error("getData: workers do not exist!")
    end
	if sum(maxProbs)<size(ProbSrcMap,1)
    	error("getMaxwellFreqParam -  maximum problems per worker are less than problems in total. ")
	end
	println("number of workers used $(length(workerList))")

	pFor  = Array(RemoteChannel,size(ProbSrcMap,1))
	probs = zeros(maximum(workerList))

	i = 1
    nextidx() = (idx=i; i+=1; idx)
	nfreq = length(freq)
	nsrc  = length(Srcs)
	# send out jobs
	@sync begin
		for p=workerList
			@async begin
				while true
					if (probs[p] >= maxProbs[p])
						break
					end
					idx = nextidx()
					if (idx > size(ProbSrcMap,1))
						break
					end
					# find src and rec on mesh
					indSrc = find(ProbSrcMap[idx,:])
					indRec = find(sum(SrcRecMap[indSrc,:],1))
					if isempty(indRec)
						error("no receiver for this src/frequency combination")
					end
					iSrc,iFreq = ind2sub((nsrc,nfreq),indSrc)
					if !all(iFreq.==iFreq[1])
						error("pFor can only handle one frequency!")
					end
					pFor[idx]  = initRemoteChannel(getMaxwellFreqParam,p,x0,n,h,Srcs[[iSrc;]],Recs[[indRec;]],freq[iFreq[1]],nf,nc,Box,linSolParam,doFV,doSE,ncells)
					wait(pFor[idx])
					probs[p] +=1
				end
			end
		end
	end
	if !(all(isdefined(pFor)))
		error("getMaxwellFreqParam - not all problems have been prepared!")
	end
	return pFor # Array of Remote Refs
end




#------ Prepare the pFor by reading Transmitter structure
function getMaxwellFreqParam(x0::Array{Float64,1},
                             n::Array{Int64,1},
                             h::Array{Float64,1},
                             TR::Union{MaxwellUtils.Transmitter, MaxwellUtils.TransmitterOmega},
                             itopo::Array,
                             meshingParam::Array,
                             linSolParam::AbstractSolver=getMUMPSsolver([],1,0),
                             doFV::Bool=true,
                             doSE=true)

    freq    = TR.omega
    
    nsmallcells   = meshingParam[1]
    mincellsize   = meshingParam[2]
    itopo         = meshingParam[3]
    depth_core    = meshingParam[4]
    mincellfactor = meshingParam[5]

    if typeof(TR) == MaxwellUtils.Transmitter
        Srcs, Recs = getSxRxFromData(TR)
    elseif typeof(TR) == MaxwellUtils.TransmitterOmega
        Srcs = TR.Srcs
        Recs = TR.Recs
    end

    if length(Srcs) != 1
        error("There should be only one source per mesh.")
    end

    # Prepare the forward modeling mesh
    M = createSmallMeshFromTX(Srcs[1], Recs, vec(h),vec(n),vec(x0),
                              nsmallcells, mincellsize, itopo,
                              depth_core, mincellfactor)


    nEx,nEy,nEz = getEdgeNumbering(M.S)

    # transmitters
    s = zeros(Complex128,sum(M.ne),length(Srcs))
    for k=1:length(Srcs)
        s[:,k] = sparse(getEdgeIntegralOfPolygonalChain(M,Srcs[k],normalize=false))
    end

    # receivers
    Obs = spzeros(Complex128,sum(M.ne),length(Recs))
    for k=1:length(Recs)
        tmp = getEdgeIntegralOfPolygonalChain(M,Recs[k],normalize=true)
        tmp = complex(tmp)
        Obs[:,k] = sparse(tmp)
        if norm(Recs[k][1,:] - Recs[k][end,:]) < 1e-16
            # for closed loops, scale by i/(omega*mu0)
            a = complex(0,1.0)/(TR.omega)/mu0
        else
            a = complex(1.0,0)
        end
        Obs[:,k] = a*Obs[:,k]
    end

    freq = TR.omega
    Sources = s

    linSolParam = getMUMPSsolver([],1,0)

	if doSE
		pFor = getMaxwellFreqParam(M,Sources,Obs,[],freq,linSolParam, sensitivityMethod=:Explicit)
	else
		pFor = getMaxwellFreqParam(M,Sources,Obs,[],freq,linSolParam)
	end

    write(STDOUT, @sprintf("    Mesh Size: %10i\n", pFor.Mesh.nc))

    return pFor
end

#---- Prepare the pFor from UBC files in parallel
function getMaxwellFreqParam(x0::Array{Float64,1},
							 n::Array{Int64,1},
							 h::Array{Float64,1},
							 TR::Union{Array{MaxwellUtils.Transmitter}, Array{MaxwellUtils.TransmitterOmega}},
							 itopo::Array,
							 meshingParam::Array,
 							 linSolParam::AbstractSolver=getMUMPSsolver([],1,0),
 							 doFV::Bool=true,
							 doSE=true)


	pFor   = Array(RemoteChannel,length(TR))
	i = 1; nextidx() = (idx=i; i+=1; idx)

	nsrc  = length(TR)
	# send out jobs
	@sync begin
		for p=workers()
			@async begin
				while true
					idx = nextidx()
					if idx > length(TR)
						break
					end
					# find src and rec on mesh
					pFor[idx]  = initRemoteChannel(getMaxwellFreqParam,
                                                                                    p,
										    vec(x0),
						   				    vec(n),
											vec(h),
											TR[idx],
											itopo,
											meshingParam,
											linSolParam,
											doFV,
											doSE)
				end
			end
		end
	end
	return pFor # Array of Remote Refs
end

#----------------------------------------------------------------

end