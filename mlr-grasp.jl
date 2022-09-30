using Random


# ==================================================================================================

mutable struct Instance
    n :: Int64
    D
    l
    avg :: Float64
    function Instance(s)
        D = zeros(1,1)
        nNodes = 0
        l = 0
        avg = 0
        if s != "STDIN"
            open(s) do f
                line = 1
                lineContent = readline(f)
                lineSplited = split(lineContent, " ")
                nNodes = parse(Int64, lineSplited[1])
                l = parse(Int64, lineSplited[2])
                D = zeros(nNodes,nNodes)
                avgRatio = 1/((nNodes^2)-nNodes)

                while ! eof(f)
                    lineContent = readline(f)
                    lineSplited = split(lineContent, " ")
                    for col = 1:nNodes
                        value = parse(Float64, lineSplited[col])
                        D[line,col] = value
                        avg += value*avgRatio
                    end

                    line += 1
                end
            end
        else
            line = 1
            lineContent = readline()
            lineSplited = split(lineContent, " ")
            nNodes = parse(Int64, lineSplited[1])
            l = parse(Int64, lineSplited[2])
            D = zeros(nNodes,nNodes)
            avgRatio = 1/((nNodes^2)-nNodes)

            while line <= nNodes
                lineContent = readline()
                lineSplited = split(lineContent, " ")
                for col = 1:nNodes
                    value = parse(Float64, lineSplited[col])
                    D[line,col] = value
                    avg += value*avgRatio
                end

                line += 1
            end
        end
        new(nNodes,D,l,avg)
    end
end

# ==================================================================================================

mutable struct Solution
    π :: Array{Int64,1} ## solução: uma permutação de cidades
    v :: Float64        ## valor da solução: distância total
    Solution(s::Solution) = new(copy(s.π),s.v)
    Solution(π,v) = new(copy(π),v)
end
Base.copy(s::Solution) = Solution(s)

# ==================================================================================================

function nearestNeighbor(I::Instance)
    s = Solution([1],0.0)
    free=trues(I.n)
    free[1]=false
    for k=2:I.n
      f=findall(free)
      j=f[findmin(I.D[s.π[end],f])[2]]
    #       s.v += I.D[s.π[end],j]
      push!(s.π,j)
      free[j]=false
    end
    for lookahead=1:I.l
        for i=1:I.n
            s.v += I.D[s.π[i],s.π[((i+lookahead-1)%(I.n))+1]]
        end
    end
#     push!(s.π,1)
    s
end

# ==================================================================================================

function randomNearestNeighbor(I::Instance,k)
    s = Solution([1],0.0)
    free=trues(I.n)
    free[1]=false
    for i=2:I.n
        f=findall(free)
        ncand=min(k,I.n-i+1)
        j=f[sortperm(I.D[s.π[end],f])[rand(1:ncand)]]
    #          s.v += I.D[s.π[end],j]
        push!(s.π,j)
        free[j]=false
    end
    for lookahead=1:I.l
        for i=1:I.n
#             println("$(s.π[i]) - $(s.π[((i+lookahead-1)%(I.n))+1])")
#             println("$(I.D[i,((i+lookahead-1)%(I.n))+1])")
            s.v += I.D[s.π[i],s.π[((i+lookahead-1)%(I.n))+1]]
        end
    end
#     push!(s.π,1)
    s
end

# ==================================================================================================

function randomNearestNeighbor(I::Instance,k)
    s = Solution([1],0.0)
    free=trues(I.n)
    free[1]=false
    for i=2:I.n
        f=findall(free)
        ncand=min(k,I.n-i+1)
        j=f[sortperm(I.D[s.π[end],f])[rand(1:ncand)]]
    #          s.v += I.D[s.π[end],j]
        push!(s.π,j)
        free[j]=false
    end
    for lookahead=1:I.l
        for i=1:I.n
#             println("$(s.π[i]) - $(s.π[((i+lookahead-1)%(I.n))+1])")
#             println("$(I.D[i,((i+lookahead-1)%(I.n))+1])")
            s.v += I.D[s.π[i],s.π[((i+lookahead-1)%(I.n))+1]]
        end
    end
#     push!(s.π,1)
    s
end

# ==================================================================================================

function calcIndex(index, n)
    result = 0
    if index <= 0
        result = n-(abs(index)%n)
    else
        result = (((index)-1)%(n))+1
    end
    result
end

function localSearchAdjacentCitiesSwap(I::Instance, s₀::Solution)
    s = copy(s₀)
    i=1
    v=0.0
    while i<I.n
        i=2
        while i<I.n
            v = s.v
            for l1=1:I.l+2
                index1 = calcIndex((i+1)-l1+1, I.n)
                for l2=1:I.l
                    index2 = calcIndex(index1+l2, I.n)
                    v -= I.D[s.π[index1],s.π[index2]]
                end
            end
            
            
            s.π[i],s.π[i+1] = s.π[i+1],s.π[i] # swap
            
            for l1=1:I.l+2
                index1 = calcIndex((i+1)-l1+1, I.n)
                for l2=1:I.l
                    index2 = calcIndex(index1+l2, I.n)
                    v += I.D[s.π[index1],s.π[index2]]
                end
            end
            
            s.π[i],s.π[i+1] = s.π[i+1],s.π[i] # swap again
            
            v<s.v && break ## first improvement
            i += 1
        end
        if i<I.n
            s.π[i],s.π[i+1] = s.π[i+1],s.π[i]
            s.v = v
        end
    end
    s
end

# ==================================================================================================

function writeInEndOfFile(content, outputFile)
    output = open(outputFile,"a")
    write(output, content);
    close(output)
end

# ==================================================================================================

# <instance input file>??
function main(args) ## args: <output file> <k> <# iterations> <seed>
#     inputFile = args[1]
    outputFile = args[1]
    k = parse(Int64, args[2])
    nIterations = parse(Int64, args[3])
    seed = parse(Int64, args[4])
    Random.seed!(seed)
    
    I = Instance("STDIN")
    
    writeInEndOfFile("Running Info | k=$(k) | #Iterations=$(nIterations) | seed=$(seed)\n", outputFile)
    
    timeNN = @elapsed begin
        b = nearestNeighbor(I)
    end
    writeInEndOfFile("Nearest Neighbor: $(b.v) | Time: $(timeNN)\n", outputFile)
    
    
    timeProgram = @elapsed begin
        writeInEndOfFile("Local Search Adjacent Cities Swap\n", outputFile)
        for i = 1:nIterations
            timeIteration = @elapsed begin
                s = randomNearestNeighbor(I,k) 
                s = localSearchAdjacentCitiesSwap(I,s)
                if s.v < b.v
                    b = copy(s)
                end
            end
            writeInEndOfFile("Iteration $i: $(b.v) | Time: $(timeIteration) \nSolution: \n$(b.π) \n", outputFile)
        end
    end
    
    
    writeInEndOfFile("Best solution: $(b.v) | Total time: $(timeProgram) \nSolution: \n$(b.π)\n\n\n\n\n\n", outputFile)
    println("Best solution: $(b.v) | Total time: $(timeProgram) \nSolution: \n$(b.π)")
end
main(ARGS)
