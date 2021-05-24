PROGRAM_VERSION = "control"
try
	using ArgParse
catch
	import Pkg
	Pkg.add("ArgParse")
	using ArgParse
end
import Statistics.median

function correction(fin, fout, X, MASK, T)
	upperx = ('a' <= X && X <= 'z') ? X - 'a' + 'A' : X
	inputText = read(fin, String)
	temp = split(inputText, ">")
	temp = temp[length.(temp) .> 0]
	temp = [split(arr, "\n") for arr in temp]
	header = [arr[1] for arr in temp]
	c = [string(arr[2:end]...) for arr in temp]
	upperc = [uppercase(str) for str in c]
	if length(c) == 0
		return
	end
	arrc = [Array{Char, 1}(str) for str in c]
	n = length(c[1])
	m = length(c)
	wo = zeros(n, m)
	# read sequences and compute per column profiles. 
	for i in 1:n
		cnt = zeros(128)
		# read a column, internally represent in upper case, and count letters
		for j in 1:m
			cnt[UInt8(upperc[j][i])] += 1
		end
		cnt[UInt8(upperx)] = 0
		cnt[UInt8('-')] = 0
		unq = length([1 for t in cnt if t > 0])
		total = sum(cnt)
		for j in 1:m
			wo[i, j] = total == 0 ? 0 : total / (unq * cnt[UInt8(upperc[j][i])]) + rand(Float64) * 0.001
		end
	end
	w1 = [wo[:,j][(arrc[j] .!= '-') .& (arrc[j] .!= X)] for j in 1:m]
	cutoff = sort([v for arr in w1 for v in arr])[end - T]
	for j in 1:m
		str = arrc[j][(arrc[j] .!= '-') .& (arrc[j] .!= X)]
		println(fout, ">" * header[j])
		strout = ""
		i = 1
		for t in 1:length(c[j])
			if c[j][t] == X || c[j][t] == '-'
				strout *= c[j][t]
			else
				strout *= (w1[j][i] > cutoff) ? MASK : str[i]
				i += 1
			end
		end
		println(fout, strout)
	end
end

function ArgParse.parse_item(::Type{Char}, x::AbstractString)
	return x[1]
end

function parse_commandline()
	s = ArgParseSettings()
	
	@add_arg_table! s begin
		"--list", "-l"
			help = "running on a list of inputs; for every three lines of the list file (ref, in, out)"
			action = :store_true
		"--mask", "-m"
			help = "the character to mask erroneous regions"
			arg_type = Char
			default = 'X'
		"--any", "-a"
			help = "the character to denote ambiguous positions or character to denote ANY in the input files"
			arg_type = Char
			default = 'X'
		"--reference", "-r"
			help = "reference with error"
			arg_type = String
			default = ""
		"input"
			help = "a fasta file as input (when -l is not set) or a list of input/output pairs (when -l is set)"
			required = true
	end
	
	return parse_args(s)
end

function main()
	println(stderr, "Version " * string(PROGRAM_VERSION))
	args = parse_commandline()

	temp = split(read(open(args["input"], "r"), String), "\n")
	temp = temp[length.(temp) .> 0]
	if args["list"] == false
		T = length([c for c in read(open(args["reference"], "r"), String) if c == args["mask"]]) - length([c for c in read(open(args["input"], "r"), String) if c == args["mask"]])
		correction(open(args["input"], "r"), stdout, args["any"], args["mask"], T)
	else
		temp = split(read(open(args["input"], "r"), String), "\n")
		temp = temp[length.(temp) .> 0]
		for i = 3:3:length(temp)
			try
				println(stderr, "Processing " * temp[i - 1] * "...")
				T = length([c for c in read(open(temp[i - 2], "r"), String) if c == args["mask"]]) - length([c for c in read(open(temp[i - 1], "r"), String) if c == args["mask"]])
				correction(open(temp[i - 1], "r"), open(temp[i], "w"), args["any"], args["mask"], T)
			catch
				println(stderr, "Error happened when processing " * temp[i - 1] * ".")
				println(stderr)
			end
		end
	end
end

main()
