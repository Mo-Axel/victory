function vech(A)
# make a vector with lower triangular matrix A (including diagonals)
A = A' # A becomes upper triangular
vech_A = A[triu(trues(size(A)))]

return vech_A

end
