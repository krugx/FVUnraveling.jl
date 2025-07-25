struct HEOMStructure
  K::Int
  T::Int
  hild::Int
  useLtrunc::Bool
  ados::Vector{Vector{Int}}        # multi-index ADOs
  idxmap::Dict{Vector{Int},Int}   # multi-index → flat index
  itvs::Vector{UnitRange{Int}}
end

struct HEOMOperator
  structure::HEOMStructure
  H::Matrix{ComplexF64}
  q::Matrix{ComplexF64}
  d::Vector{ComplexF64}
  gamma::Vector{ComplexF64}

end

struct CachedOps
  Id_op::SparseMatrixCSC
  L_op::SparseMatrixCSC
  Q_op::SparseMatrixCSC
  Q_L_op::SparseMatrixCSC
  Q_R_op::SparseMatrixCSC
  function CachedOps(H, q)
    id = sparse(diagm(ones(size(H, 1))))
    Id_op = kron(id', id)
    L_op = kron(id', H) - kron(H', id)
    Id_op = kron(id', id)
    Q_op = kron(id', q) - kron(q', id)
    Q_L_op = kron(id', q)
    Q_R_op = kron(q', id)

    new(Id_op, L_op, Q_op, Q_L_op, Q_R_op)
  end
end

struct HEOMPropagator
  mat::SparseMatrixCSC
  structure::HEOMStructure
  function HEOMPropagator(op::HEOMOperator)
    new(build_matrix(op), op.structure)
  end
  function HEOMPropagator(structure, H, q, d, gamma)
    op = HEOMOperator(structure, H, q, d, gamma)
    new(build_matrix(op), structure)
  end
end

function build_heom_structure(K::Int, T::Int, hild::Int; useLtrunc::Bool=true)
  ados = Vector{Vector{Int}}()
  idxmap = Dict{Vector{Int},Int}()
  itvs = Vector{UnitRange{Int}}()

  counter = 0
  current = zeros(Int, 2K)

  function traverse(i::Int, sum::Int)
    if i > 2K
      r = copy(current)
      push!(ados, r)
      idxmap[r] = counter
      push!(itvs, counter*hild^2+1:(counter+1)*hild^2)
      counter += 1
      return
    end
    for val in 0:T-1
      if useLtrunc && (sum + val > T - 1)
        break
      end
      current[i] = val
      traverse(i + 1, sum + val)
    end
  end
  traverse(1, 0)

  return HEOMStructure(K, T, hild, useLtrunc, ados, idxmap, itvs)
end

function diag_ado(r, op::HEOMOperator, cops::CachedOps)
  st = op.structure

  ado_op = -im * cops.L_op
  ado_op += -sum(r[1:st.K] .* op.gamma[1:st.K] + r[st.K+1:2st.K] .* conj.(op.gamma)) * cops.Id_op

  return ado_op
end

function shift_ado(r, j::Int, offset::Int, direction::Int, op::HEOMOperator, cops::CachedOps)
  st = op.structure

  d_ext = vcat(op.d, conj.(op.d))

  if direction == -1
    return -im * sqrt(r[j+offset] * d_ext[j+offset]) * cops.Q_op
  elseif direction == +1
    if offset == 0
      return -im * sqrt((r[j+offset] + 1) * d_ext[j+offset]) * cops.Q_L_op
    elseif offset == st.K
      return +im * sqrt((r[j+offset] + 1) * d_ext[j+offset]) * cops.Q_R_op
    end
  end
end

function build_matrix(op::HEOMOperator)::SparseMatrixCSC
  st = op.structure

  offset_lst = [0, st.K]
  direction_lst = [-1, +1]

  dim = length(st.ados) * st.hild^2
  heom_matrix = spzeros(ComplexF64, dim, dim)

  cops = CachedOps(op.H, op.q)

  for i in eachindex(st.ados)
    r = st.ados[i]
    itv = st.itvs[i]

    ## Diag elements
    heom_matrix[itv, itv] = diag_ado(r, op, cops)

    ## Offdiagonal elements
    for j in 1:st.K
      for offset in offset_lst
        for direction in direction_lst
          r_shift = copy(r)
          r_shift[offset+j] += direction
          if haskey(st.idxmap, r_shift)
            itv_shift = st.itvs[st.idxmap[r_shift]+1]
            heom_matrix[itv_shift, itv] += shift_ado(r, j, offset, direction, op, cops)
          end
        end
      end
    end
  end

  return heom_matrix
end

function check_stability(prop::HEOMPropagator)
  st = prop.structure
  dim = length(st.ados) * st.hild^2

  v = zeros(ComplexF64, dim)
  # initial guess for steady state in convergence
  v[1:st.hild^2] = reshape(1 / 2 * [1 1-im; 1+im 1], (st.hild^2, 1))
  lambda, v_lambda = eigs(prop.mat, nev=1; which=:LR, maxiter=2000, check=1, v0=v)

  lambda = real(lambda[1])
  v_lambda = sparsevec(v_lambda)

  ## Check if stable
  println("max({Re(lambda_i)}) = $lambda \n")
  if lambda <= 1e-6
    rho_SS = reshape(v_lambda[1:st.hild^2], st.hild, st.hild)
    phys_SS = rho_SS / tr(rho_SS)
    ado_SS = v_lambda / tr(rho_SS)

    println("Physical steady state rho=")
    display(rho_SS)
    println("\n")
    println("Steady state of ADO state = ")
    display(v_lambda)
    return phys_SS, ado_SS
  else
    println("Instable HEOM dynamics!")
    return -1
  end
end
